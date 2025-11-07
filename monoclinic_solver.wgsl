// monoclinic_solver.wgsl

// === Structs ===
struct RawMonoSolution {
    a: f32,
    b: f32,
    c: f32,
    beta: f32,
}

// === Type Aliases ===
alias Vec4 = vec4<f32>;
alias Mat4x4 = mat4x4<f32>; 

// === Bindings ===
@group(0) @binding(0) var<storage, read> q_obs: array<f32>;
@group(0) @binding(1) var<storage, read> hkl_basis: array<f32>; // [h,k,l,pad]
@group(0) @binding(2) var<storage, read> peak_combos: array<u32>; // [i,j,k,l]
@group(0) @binding(3) var<storage, read> hkl_combos: array<u32>; // [n1,n2,n3,n4]
@group(0) @binding(4) var<storage, read_write> solution_counter: atomic<u32>;
@group(0) @binding(5) var<storage, read_write> results_list: array<RawMonoSolution>;

// *** NEW ***
struct Config { z_offset: u32, };
@group(0) @binding(6) var<uniform> config: Config;
// *** END NEW ***

// === Constants ===
const PI: f32 = 3.1415926535;
const RAD: f32 = PI / 180.0;
const DEG: f32 = 180.0 / PI;
const WORKGROUP_SIZE_Y: u32 = 8u; // Must match @workgroup_size
const MAX_Y_WORKGROUPS: u32 = 65535u; // WebGPU limit

// ... (Helper Functions solve4x4, extractCell, plausible all stay the same) ...

// === Helper Functions ===

// 4x4 Matrix Solve (Gaussian Elimination)
fn solve4x4(A: Mat4x4, b: Vec4) -> Vec4 {
    var M: Mat4x4 = A; 
    var v: Vec4 = b;
    let n: u32 = 4u;

    for (var i: u32 = 0u; i < n; i = i + 1u) {
        let pivot: f32 = M[i][i]; 
        if (abs(pivot) < 1e-10) { 
            return Vec4(0.0, 0.0, 0.0, 0.0); 
        }

        for (var r: u32 = i + 1u; r < n; r = r + 1u) {
            let fac: f32 = M[r][i] / pivot;
            M[r] = M[r] - (fac * M[i]);
            v[r] = v[r] - (fac * v[i]);
        }
    }

    var x: Vec4;
    for (var i_s: i32 = 3; i_s >= 0; i_s = i_s - 1) {
        let i: u32 = u32(i_s);
        var s: f32 = v[i];
        for (var j: u32 = i + 1u; j < n; j = j + 1u) {
            s = s - M[i][j] * x[j];
        }
        x[i] = s / M[i][i];
    }
    return x;
}

// Extract cell from 4 params
fn extractCell(params: Vec4) -> RawMonoSolution {
    let A: f32 = params[0]; // h*h
    let B: f32 = params[1]; // k*k
    let C: f32 = params[2]; // l*l
    let D: f32 = params[3]; // h*l

    if (A <= 0.0 || B <= 0.0 || C <= 0.0 || (A != A) || (B != B) || (C != C)) {
        return RawMonoSolution();
    }
    if (D*D >= 4.0 * A * C) {
        return RawMonoSolution();
    }
    let cosBeta_calc = -D / (2.0 * sqrt(A*C));
    if (abs(cosBeta_calc) >= 1.0) {
        return RawMonoSolution();
    }
    var beta_calc = acos(cosBeta_calc) * DEG;
    if (beta_calc < 90.0) {
        beta_calc = 180.0 - beta_calc;
    }
    if (beta_calc < 90.0 || beta_calc > 150.0) {
        return RawMonoSolution();
    }
    let sinBetaSq = sin(beta_calc * RAD) * sin(beta_calc * RAD);
    if (sinBetaSq <= 1e-6) {
        return RawMonoSolution();
    }
    return RawMonoSolution(
        1.0 / sqrt(A * sinBetaSq), // a
        1.0 / sqrt(B),             // b
        1.0 / sqrt(C * sinBetaSq), // c
        beta_calc                  // beta
    );
}

// Simple plausibility check
fn plausible(cell: RawMonoSolution) -> bool {
    if (cell.a < 2.0 || cell.a > 50.0 || (cell.a != cell.a)) { return false; }
    if (cell.b < 2.0 || cell.b > 50.0 || (cell.b != cell.b)) { return false; }
    if (cell.c < 2.0 || cell.c > 50.0 || (cell.c != cell.c)) { return false; }
    return true;
}

// === Main Kernel ===
@compute @workgroup_size(8, WORKGROUP_SIZE_Y, 1)
fn main(
    @builtin(global_invocation_id) global_id: vec3<u32>,
    @builtin(num_workgroups) num_workgroups: vec3<u32>
) { 
    
    // 1. Get indices for this thread
    let peak_combo_idx: u32 = global_id.x;

    // *** MODIFIED ***
    // Calculate HKL index based on Z-offset from uniforms
    let hkl_threads_per_z_slice = MAX_Y_WORKGROUPS * WORKGROUP_SIZE_Y;
    let z_index = global_id.z + config.z_offset; // global_id.z will be 0
    let hkl_combo_idx: u32 = z_index * hkl_threads_per_z_slice + global_id.y;
    // *** END MODIFIED ***


    // Bounds check
    let num_peak_combos = arrayLength(&peak_combos) / 4u;
    let num_hkl_combos = arrayLength(&hkl_combos) / 4u;

    if (peak_combo_idx >= num_peak_combos || 
        hkl_combo_idx >= num_hkl_combos) {
        return;
    }

    // 2. Get q_obs vector (NO PERMUTATION)
    let p_offset = peak_combo_idx * 4u;
    let q_vec = Vec4(
        q_obs[peak_combos[p_offset + 0u]],
        q_obs[peak_combos[p_offset + 1u]],
        q_obs[peak_combos[p_offset + 2u]],
        q_obs[peak_combos[p_offset + 3u]]
    );

    // 3. Get HKLs and build M matrix
    let h_offset = hkl_combo_idx * 4u;
    var M_cols: array<vec4<f32>, 4>;
    for (var i: u32 = 0u; i < 4u; i = i + 1u) {
        let hkl_idx = hkl_combos[h_offset + i];
        let h = hkl_basis[hkl_idx * 4u + 0u];
        let k = hkl_basis[hkl_idx * 4u + 1u];
        let l = hkl_basis[hkl_idx * 4u + 2u];
        
        M_cols[0][i] = h*h;
        M_cols[1][i] = k*k;
        M_cols[2][i] = l*l;
        M_cols[3][i] = h*l;
    }
    var M = mat4x4<f32>(M_cols[0], M_cols[1], M_cols[2], M_cols[3]);
    M = transpose(M); 

    // 4. Solve the system
    let fit_params = solve4x4(M, q_vec);
    if (abs(fit_params[0]) < 1e-10) { return; } // Solve failed

    // 5. Extract cell and check it
    let cell = extractCell(fit_params);
    if (!plausible(cell)) {
        return;
    }

    // 6. Save the plausible solution
    let solution_index = atomicAdd(&solution_counter, 1u);
    if (solution_index < 50000u) { // maxSolutions
        results_list[solution_index] = cell;
    }
}