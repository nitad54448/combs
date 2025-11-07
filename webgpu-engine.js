// webgpu-engine.js
class WebGPUEngine {
    constructor() {
        this.device = null;
        this.adapter = null; 
        this.shaderModule = null;
        this.pipeline = null;
    }

    // 1. Initialize WebGPU
    async init() {
        if (!navigator.gpu) {
            throw new Error("WebGPU not supported on this browser.");
        }
        this.adapter = await navigator.gpu.requestAdapter(); 
        if (!this.adapter) { 
            throw new Error("No compatible GPUAdapter found.");
        }
        this.device = await this.adapter.requestDevice(); 
        return true;
    }

    // 2. Load the WGSL Shader
    async loadShader(url) {
        const response = await fetch(url);
        const shaderCode = await response.text();
        this.shaderModule = this.device.createShaderModule({ code: shaderCode });
    }

    // 3. Create the compute pipeline
    createPipeline(entryPoint = "main") {
        this.pipeline = this.device.createComputePipeline({
            layout: 'auto',
            compute: {
                module: this.shaderModule,
                entryPoint: entryPoint,
            },
        });
    }

    // 4. Helper to create a buffer and write data to it
    createBuffer(data, usage) {
        const buffer = this.device.createBuffer({
            size: data.byteLength,
            usage: usage | GPUBufferUsage.COPY_DST,
            mappedAtCreation: true,
        });
        new data.constructor(buffer.getMappedRange()).set(data);
        buffer.unmap();
        return buffer;
    }

    // 5. Helper to create a buffer for reading results back
    createReadBuffer(size) {
        return this.device.createBuffer({
            size: size,
            usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST,
        });
    }

    // 6. Helper to create a buffer for GPU-only storage
    createStorageBuffer(size) {
        return this.device.createBuffer({
            size: size,
            usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST | GPUBufferUsage.COPY_SRC,
        });
    }

    // 7. This is the main function you'll call
    async runTriclinicSolver(
        qObsArray,        // Float32Array of q_obs values
        hklBasisArray,    // Float32Array of [h,k,l,pad, h,k,l,pad, ...]
        peakCombos,       // Uint32Array of [i,j,k,l,m,n, i,j,k,l,m,n, ...]
        hklCombos,        // Uint32Array of [n1,n2,n3,n4,n5,n6, ...]
        progressCallback  // *** NEW ***
    ) {
        if (!this.pipeline) {
            throw new Error("Pipeline not created. Call createPipeline() first.");
        }
        if (!this.adapter) { 
            throw new Error("Engine not initialized. Call init() first.");
        }

        // --- Create Buffers ---
        const maxSolutions = 50000; 
        const solutionStructSize = 6 * 4; 
        
        const qObsBuffer = this.createBuffer(qObsArray, GPUBufferUsage.STORAGE);
        const hklBasisBuffer = this.createBuffer(hklBasisArray, GPUBufferUsage.STORAGE);
        const peakCombosBuffer = this.createBuffer(peakCombos, GPUBufferUsage.STORAGE);
        const hklCombosBuffer = this.createBuffer(hklCombos, GPUBufferUsage.STORAGE);

        const counterBuffer = this.createStorageBuffer(4); 
        const resultsBuffer = this.createStorageBuffer(maxSolutions * solutionStructSize);
        const counterReadBuffer = this.createReadBuffer(4);
        const resultsReadBuffer = this.createReadBuffer(maxSolutions * solutionStructSize);

        // *** NEW *** Uniform buffer for z_offset
        const configBuffer = this.device.createBuffer({
            size: 4, // one u32
            usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
        });
        const configUniformArray = new Uint32Array(1); // JS-side array to write from
        // *** END NEW ***

        // --- Create Bind Group (once) ---
        const bindGroup = this.device.createBindGroup({
            layout: this.pipeline.getBindGroupLayout(0),
            entries: [
                { binding: 0, resource: { buffer: qObsBuffer } },
                { binding: 1, resource: { buffer: hklBasisBuffer } },
                { binding: 2, resource: { buffer: peakCombosBuffer } },
                { binding: 3, resource: { buffer: hklCombosBuffer } },
                { binding: 4, resource: { buffer: counterBuffer } }, 
                { binding: 5, resource: { buffer: resultsBuffer } },
                { binding: 6, resource: { buffer: configBuffer } }, // *** NEW ***
            ],
        });

        // --- *** NEW Chunked Dispatch Logic *** ---
        
        const numPeakCombos = peakCombos.length / 6;
        const numHklCombos = hklCombos.length / 6;
        const workgroupSizeX = 4; // From WGSL @workgroup_size
        const workgroupSizeY = 4; // From WGSL @workgroup_size

        const workgroupsX = Math.ceil(numPeakCombos / workgroupSizeX);
        const totalHklWorkgroups = Math.ceil(numHklCombos / workgroupSizeY);
        const maxDimY = this.adapter.limits.maxComputeWorkgroupsPerDimension || 65535;

        // We will dispatch in Z-chunks, where each chunk has a full Y-dimension
        const workgroupsY = maxDimY; 
        const totalWorkgroupsZ = Math.ceil(totalHklWorkgroups / workgroupsY);

        for (let z_chunk = 0; z_chunk < totalWorkgroupsZ; z_chunk++) {
            // 1. Update uniform buffer with the current Z-offset
            configUniformArray[0] = z_chunk;
            this.device.queue.writeBuffer(configBuffer, 0, configUniformArray);

            // 2. Calculate workgroups for *this* chunk
            // If this is the last chunk, it might not be a full 'maxDimY'
            const hklWorkgroupsInThisChunk = (z_chunk === totalWorkgroupsZ - 1) 
                ? (totalHklWorkgroups % workgroupsY || workgroupsY) 
                : workgroupsY;

            // 3. Create and submit commands for this chunk
            const commandEncoder = this.device.createCommandEncoder();
            const passEncoder = commandEncoder.beginComputePass();
            passEncoder.setPipeline(this.pipeline);
            passEncoder.setBindGroup(0, bindGroup);
            // Dispatch (X, Y_chunk_size, 1)
            passEncoder.dispatchWorkgroups(workgroupsX, hklWorkgroupsInThisChunk, 1); 
            passEncoder.end();
            
            // We also poll the counter *after* this chunk
            commandEncoder.copyBufferToBuffer(counterBuffer, 0, counterReadBuffer, 0, 4);
            
            this.device.queue.submit([commandEncoder.finish()]);
            await this.device.queue.onSubmittedWorkDone(); // Wait for this chunk

            // 4. Report Progress
            if (progressCallback) {
                await counterReadBuffer.mapAsync(GPUMapMode.READ);
                const numSolutions = new Uint32Array(counterReadBuffer.getMappedRange())[0];
                counterReadBuffer.unmap();
                // Send progress (0.0 to 1.0) and current solution count
                progressCallback((z_chunk + 1) / totalWorkgroupsZ, numSolutions);
            }
        }
        // --- *** End of New Logic *** ---

        // After all chunks are done, copy the *final* results
        const finalEncoder = this.device.createCommandEncoder();
        finalEncoder.copyBufferToBuffer(counterBuffer, 0, counterReadBuffer, 0, 4);
        finalEncoder.copyBufferToBuffer(resultsBuffer, 0, resultsReadBuffer, 0, resultsBuffer.size);
        this.device.queue.submit([finalEncoder.finish()]);
        await this.device.queue.onSubmittedWorkDone();


        // --- Read Results ---
        await counterReadBuffer.mapAsync(GPUMapMode.READ);
        const numSolutions = new Uint32Array(counterReadBuffer.getMappedRange())[0];
        counterReadBuffer.unmap();

        const potentialCells = [];
        if (numSolutions > 0) {
            await resultsReadBuffer.mapAsync(GPUMapMode.READ);
            const rawResults = new Float32Array(resultsReadBuffer.getMappedRange());
            
            for (let i = 0; i < Math.min(numSolutions, maxSolutions); i++) {
                const offset = i * 6;
                potentialCells.push({
                    a: rawResults[offset + 0],
                    b: rawResults[offset + 1],
                    c: rawResults[offset + 2],
                    alpha: rawResults[offset + 3],
                    beta: rawResults[offset + 4],
                    gamma: rawResults[offset + 5],
                    system: 'triclinic',
                });
            }
            resultsReadBuffer.unmap();
        }

        // --- Cleanup ---
        qObsBuffer.destroy();
        hklBasisBuffer.destroy();
        peakCombosBuffer.destroy();
        hklCombosBuffer.destroy();
        counterBuffer.destroy();
        resultsBuffer.destroy();
        counterReadBuffer.destroy();
        resultsReadBuffer.destroy();
        configBuffer.destroy(); // *** NEW ***

        return potentialCells;
    }

    // 8. This is the main function you'll call
    async runMonoclinicSolver(
        qObsArray,        // Float32Array of q_obs values
        hklBasisArray,    // Float32Array of [h,k,l,pad, h,k,l,pad, ...]
        peakCombos,       // Uint32Array of [i,j,k,l, i,j,k,l, ...]
        hklCombos,        // Uint32Array of [n1,n2,n3,n4, ...]
        progressCallback  // *** NEW ***
    ) {
        if (!this.pipeline) {
            throw new Error("Pipeline not created. Call createPipeline() first.");
        }
        if (!this.adapter) { 
            throw new Error("Engine not initialized. Call init() first.");
        }

        // --- Create Buffers ---
        const maxSolutions = 50000;
        const solutionStructSize = 4 * 4; 
        
        const qObsBuffer = this.createBuffer(qObsArray, GPUBufferUsage.STORAGE);
        const hklBasisBuffer = this.createBuffer(hklBasisArray, GPUBufferUsage.STORAGE);
        const peakCombosBuffer = this.createBuffer(peakCombos, GPUBufferUsage.STORAGE);
        const hklCombosBuffer = this.createBuffer(hklCombos, GPUBufferUsage.STORAGE);

        const counterBuffer = this.createStorageBuffer(4);
        const resultsBuffer = this.createStorageBuffer(maxSolutions * solutionStructSize);
        const counterReadBuffer = this.createReadBuffer(4);
        const resultsReadBuffer = this.createReadBuffer(maxSolutions * solutionStructSize);

        // *** NEW *** Uniform buffer for z_offset
        const configBuffer = this.device.createBuffer({
            size: 4, // one u32
            usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST,
        });
        const configUniformArray = new Uint32Array(1); // JS-side array to write from
        // *** END NEW ***

        // --- Create Bind Group (once) ---
        const bindGroup = this.device.createBindGroup({
            layout: this.pipeline.getBindGroupLayout(0),
            entries: [
                { binding: 0, resource: { buffer: qObsBuffer } },
                { binding: 1, resource: { buffer: hklBasisBuffer } },
                { binding: 2, resource: { buffer: peakCombosBuffer } },
                { binding: 3, resource: { buffer: hklCombosBuffer } },
                { binding: 4, resource: { buffer: counterBuffer } },
                { binding: 5, resource: { buffer: resultsBuffer } },
                { binding: 6, resource: { buffer: configBuffer } }, // *** NEW ***
            ],
        });

        // --- *** NEW Chunked Dispatch Logic *** ---
        
        const numPeakCombos = peakCombos.length / 4;
        const numHklCombos = hklCombos.length / 4;
        const workgroupSizeX = 8; // From WGSL @workgroup_size
        const workgroupSizeY = 8; // From WGSL @workgroup_size

        const workgroupsX = Math.ceil(numPeakCombos / workgroupSizeX);
        const totalHklWorkgroups = Math.ceil(numHklCombos / workgroupSizeY);
        const maxDimY = this.adapter.limits.maxComputeWorkgroupsPerDimension || 65535;

        const workgroupsY = maxDimY; 
        const totalWorkgroupsZ = Math.ceil(totalHklWorkgroups / workgroupsY);

        for (let z_chunk = 0; z_chunk < totalWorkgroupsZ; z_chunk++) {
            // 1. Update uniform buffer
            configUniformArray[0] = z_chunk;
            this.device.queue.writeBuffer(configBuffer, 0, configUniformArray);

            // 2. Calculate workgroups for this chunk
            const hklWorkgroupsInThisChunk = (z_chunk === totalWorkgroupsZ - 1) 
                ? (totalHklWorkgroups % workgroupsY || workgroupsY) 
                : workgroupsY;

            // 3. Create and submit commands
            const commandEncoder = this.device.createCommandEncoder();
            const passEncoder = commandEncoder.beginComputePass();
            passEncoder.setPipeline(this.pipeline);
            passEncoder.setBindGroup(0, bindGroup);
            passEncoder.dispatchWorkgroups(workgroupsX, hklWorkgroupsInThisChunk, 1);
            passEncoder.end();
            
            commandEncoder.copyBufferToBuffer(counterBuffer, 0, counterReadBuffer, 0, 4);
            
            this.device.queue.submit([commandEncoder.finish()]);
            await this.device.queue.onSubmittedWorkDone();

            // 4. Report Progress
            if (progressCallback) {
                await counterReadBuffer.mapAsync(GPUMapMode.READ);
                const numSolutions = new Uint32Array(counterReadBuffer.getMappedRange())[0];
                counterReadBuffer.unmap();
                progressCallback((z_chunk + 1) / totalWorkgroupsZ, numSolutions);
            }
        }
        // --- *** End of New Logic *** ---

        // After all chunks are done, copy the *final* results
        const finalEncoder = this.device.createCommandEncoder();
        finalEncoder.copyBufferToBuffer(counterBuffer, 0, counterReadBuffer, 0, 4);
        finalEncoder.copyBufferToBuffer(resultsBuffer, 0, resultsReadBuffer, 0, resultsBuffer.size);
        this.device.queue.submit([finalEncoder.finish()]);
        await this.device.queue.onSubmittedWorkDone();


        // --- Read Results ---
        await counterReadBuffer.mapAsync(GPUMapMode.READ);
        const numSolutions = new Uint32Array(counterReadBuffer.getMappedRange())[0];
        counterReadBuffer.unmap();

        const potentialCells = [];
        if (numSolutions > 0) {
            await resultsReadBuffer.mapAsync(GPUMapMode.READ);
            const rawResults = new Float32Array(resultsReadBuffer.getMappedRange());
            
            for (let i = 0; i < Math.min(numSolutions, maxSolutions); i++) {
                const offset = i * 4;
                potentialCells.push({
                    a: rawResults[offset + 0],
                    b: rawResults[offset + 1],
                    c: rawResults[offset + 2],
                    beta: rawResults[offset + 3],
                    system: 'monoclinic',
                });
            }
            resultsReadBuffer.unmap();
        }

        // --- Cleanup ---
        qObsBuffer.destroy();
        hklBasisBuffer.destroy();
        peakCombosBuffer.destroy();
        hklCombosBuffer.destroy();
        counterBuffer.destroy();
        resultsBuffer.destroy();
        counterReadBuffer.destroy();
        resultsReadBuffer.destroy();
        configBuffer.destroy(); // *** NEW ***

        return potentialCells;
    }
}