// triclinic_solver.wgsl

// === Structs ===
struct RawSolution {
    a: f32,
    b: f32,
    c: f32,
    alpha: f32,
    beta: f32,
    gamma: f32,
}

// === Type Aliases ===
alias Vec6 = array<f32, 6>;
alias Mat6x6 = array<f32, 36>; // A flat 6x6 matrix, stored row-major

// === Bindings ===
@group(0) @binding(0) var<storage, read> q_obs: array<f32>;
@group(0) @binding(1) var<storage, read> hkl_basis: array<f32>; // [h,k,l,pad]
@group(0) @binding(2) var<storage, read> peak_combos: array<u32>; // [i,j,k,l,m,n]
@group(0) @binding(3) var<storage, read> hkl_combos: array<u32>; // [n1,n2,n3,n4,n5,n6]
@group(0) @binding(4) var<storage, read_write> solution_counter: atomic<u32>;
@group(0) @binding(5) var<storage, read_write> results_list: array<RawSolution>;

// *** NEW ***
struct Config { z_offset: u32, };
@group(0) @binding(6) var<uniform> config: Config;
// *** END NEW ***

// === Constants ===
const PI: f32 = 3.1415926535;
const RAD: f32 = PI / 180.0;
const DEG: f32 = 180.0 / PI;
const WORKGROUP_SIZE_Y: u32 = 4u; // Must match @workgroup_size
const MAX_Y_WORKGROUPS: u32 = 65535u; // WebGPU limit

// ... (Helper Functions solve6x6, invert3x3, extractCell, plausible all stay the same) ...

// === Helper Functions (Ported from your JS) ===

// 6x6 Matrix Solve (Gaussian Elimination)
fn solve6x6(A: Mat6x6, b: Vec6) -> Vec6 {
    var M: Mat6x6 = A;
    var v: Vec6 = b;
    let n: u32 = 6u;

    for (var i: u32 = 0u; i < n; i = i + 1u) {
        let pivot: f32 = M[i * n + i]; // M[i][i]
        if (abs(pivot) < 1e-10) { 
            return Vec6(0.0, 0.0, 0.0, 0.0, 0.0, 0.0); 
        }

        for (var r: u32 = i + 1u; r < n; r = r + 1u) {
            let fac: f32 = M[r * n + i] / pivot; // M[r][i] / pivot
            for (var c: u32 = i; c < n; c = c + 1u) {
                M[r * n + c] = M[r * n + c] - fac * M[i * n + c]; // M[r][c] -= fac * M[i][c]
            }
            v[r] = v[r] - fac * v[i];
        }
    }

    var x: Vec6;
    for (var i_s: i32 = 5; i_s >= 0; i_s = i_s - 1) {
        let i: u32 = u32(i_s);
        var s: f32 = v[i];
        for (var j: u32 = i + 1u; j < n; j = j + 1u) {
            s = s - M[i * n + j] * x[j]; // s -= M[i][j] * x[j]
        }
        x[i] = s / M[i * n + i]; // x[i] = s / M[i][i]
    }
    return x;
}

// Invert 3x3
fn invert3x3(M: mat3x3<f32>) -> mat3x3<f32> {
    let det = determinant(M);
    if (abs(det) < 1e-14) { return mat3x3<f32>(); }

    let invDet: f32 = 1.0 / det;
    var inv: mat3x3<f32>;
    
    inv[0][0] = (M[1][1] * M[2][2] - M[1][2] * M[2][1]) * invDet;
    inv[0][1] = (M[0][2] * M[2][1] - M[0][1] * M[2][2]) * invDet;
    inv[0][2] = (M[0][1] * M[1][2] - M[0][2] * M[1][1]) * invDet;
    inv[1][0] = (M[1][2] * M[2][0] - M[1][0] * M[2][2]) * invDet;
    inv[1][1] = (M[0][0] * M[2][2] - M[0][2] * M[2][0]) * invDet;
    inv[1][2] = (M[0][2] * M[1][0] - M[0][0] * M[1][2]) * invDet;
    inv[2][0] = (M[1][0] * M[2][1] - M[1][1] * M[2][0]) * invDet;
    inv[2][1] = (M[0][1] * M[2][0] - M[0][0] * M[2][1]) * invDet;
    inv[2][2] = (M[0][0] * M[1][1] - M[0][1] * M[1][0]) * invDet;

    return inv;
}

// Extract cell from 6 params (G* -> G -> cell)
fn extractCell(params: Vec6) -> RawSolution {
    let p1: f32 = params[0];
    let p2: f32 = params[1];
    let p3: f32 = params[2];
    let p4: f32 = params[3];
    let p5: f32 = params[4];
    let p6: f32 = params[5];

    let G_star = mat3x3<f32>(
        vec3<f32>(p1, p6/2.0, p5/2.0),
        vec3<f32>(p6/2.0, p2, p4/2.0),
        vec3<f32>(p5/2.0, p4/2.0, p3)
    );
    
    let G = invert3x3(G_star);
    if (abs(G[0][0]) < 1e-10) { return RawSolution(); } // Invalid

    let a: f32 = sqrt(G[0][0]);
    let b: f32 = sqrt(G[1][1]);
    let c: f32 = sqrt(G[2][2]);

    if (a < 1e-6 || b < 1e-6 || c < 1e-6 || (a != a) || (b != b) || (c != c)) {
        return RawSolution(); // Invalid
    }

    let alpha_cos = clamp(G[1][2] / (b*c), -1.0, 1.0);
    let beta_cos = clamp(G[0][2] / (a*c), -1.0, 1.0);
    let gamma_cos = clamp(G[0][1] / (a*b), -1.0, 1.0);

    return RawSolution(
        a, b, c,
        acos(alpha_cos) * DEG,
        acos(beta_cos) * DEG,
        acos(gamma_cos) * DEG
    );
}

// Simple plausibility check
fn plausible(cell: RawSolution) -> bool {
    if (cell.a < 2.0 || cell.a > 50.0 || (cell.a != cell.a)) { return false; }
    if (cell.b < 2.0 || cell.b > 50.0 || (cell.b != cell.b)) { return false; }
    if (cell.c < 2.0 || cell.c > 50.0 || (cell.c != cell.c)) { return false; }
    if (cell.alpha < 60.0 || cell.alpha > 150.0 || (cell.alpha != cell.alpha)) { return false; }
    if (cell.beta < 60.0 || cell.beta > 150.0 || (cell.beta != cell.beta)) { return false; }
    if (cell.gamma < 60.0 || cell.gamma > 150.0 || (cell.gamma != cell.gamma)) { return false; }
    return true;
}

// === Main Kernel ===
@compute @workgroup_size(4, WORKGROUP_SIZE_Y, 1)
fn main(
    @builtin(global_invocation_id) global_id: vec3<u32>,
    @builtin(num_workgroups) num_workgroups: vec3<u32>
) {
    
    // 1. Get indices for this thread
    let peak_combo_idx: u32 = global_id.x;
    
    // *** MODIFIED ***
    // Calculate HKL index based on Z-offset from uniforms
    // We hardcode MAX_Y_WORKGROUPS to match the JS dispatch logic
    let hkl_threads_per_z_slice = MAX_Y_WORKGROUPS * WORKGROUP_SIZE_Y;
    // global_id.z will be 0 (since we dispatch Z=1), config.z_offset is our loop index
    let z_index = global_id.z + config.z_offset; 
    let hkl_combo_idx: u32 = z_index * hkl_threads_per_z_slice + global_id.y;
    // *** END MODIFIED ***

    // Bounds check
    let num_peak_combos = arrayLength(&peak_combos) / 6u;
    let num_hkl_combos = arrayLength(&hkl_combos) / 6u;

    if (peak_combo_idx >= num_peak_combos || 
        hkl_combo_idx >= num_hkl_combos) {
        return;
    }

    // 2. Get q_obs vector
    let p_offset = peak_combo_idx * 6u;
    let q_vec = Vec6(
        q_obs[peak_combos[p_offset + 0u]],
        q_obs[peak_combos[p_offset + 1u]],
        q_obs[peak_combos[p_offset + 2u]],
        q_obs[peak_combos[p_offset + 3u]],
        q_obs[peak_combos[p_offset + 4u]],
        q_obs[peak_combos[p_offset + 5u]]
    );

    // 3. Get HKLs and build M matrix
    let h_offset = hkl_combo_idx * 6u;
    var M: Mat6x6; // This is our flat array<f32, 36>
    
    for (var i: u32 = 0u; i < 6u; i = i + 1u) {
        let hkl_idx = hkl_combos[h_offset + i];
        let h = hkl_basis[hkl_idx * 4u + 0u];
        let k = hkl_basis[hkl_idx * 4u + 1u];
        let l = hkl_basis[hkl_idx * 4u + 2u];
        
        let row_offset = i * 6u;
        M[row_offset + 0u] = h*h;
        M[row_offset + 1u] = k*k;
        M[row_offset + 2u] = l*l;
        M[row_offset + 3u] = k*l;
        M[row_offset + 4u] = h*l;
        M[row_offset + 5u] = h*k;
    }

    // 4. Solve the system
    let fit_params = solve6x6(M, q_vec);
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