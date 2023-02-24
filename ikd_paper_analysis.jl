using Revise
using QSM: ikd, rts, nltv, tikh, tsvd, tkd, ndi
using QSM: padarray!, padfastfft, dipole_kernel

using MAT: matread
using NIfTI: niread, niwrite, NIVolume
using DelimitedFiles: writedlm

using Images: dilate, erode, mse
using ImageQualityIndexes: assess, PSNR, XSIM
using LinearAlgebra: norm

using FFTW: plan_fft
using LinearOperators: LinearOperator
using SparsityOperators: linearOperator # for Wavelet (db2)
using RegularizedLeastSquares: Regularization, createLinearSolver, solve

# constants
γ = 267.52 # gyromagnetic ratio
dataset = "Challenge"
sweepBandLimitFlag = false
sweepMaskFlag = false
tuningFlag = false

if dataset == "Challenge"
    mag = niread("data/ChallengeSim2SNR1/Magnitude.nii.gz")
    ϕ = niread("data/ChallengeSim2SNR1/Frequency.nii.gz").raw # Hz
    mask = Bool.(niread("data/ChallengeSim2SNR1/Brain.nii.gz").raw)
    χ = niread("data/ChallengeSim2SNR1/Chi.nii.gz").raw

    params = matread("data/ChallengeSim2SNR1/SimulationParameters.mat")

    B0 = params["SimParams"]["B0"]
    vsz = tuple(params["SimParams"]["Res"]...)
    TEs = params["SeqParams"]["TE"]
    sz = size(ϕ)

    ϕ .*= TEs[1]
    # ϕ .*= 2*π # hz => rad
    # ϕ .*= 1e-3 # rad => ppm
    if !(isdir("./results/$(dataset)"))
        mkdir("./results/$(dataset)")
    end
elseif dataset == "AnitaBrain5"
    mag = niread("data/AnitaCoverageBrain5/Thickness_1mm_magnitude.nii.gz")
    ϕ = niread("data/AnitaCoverageBrain5/Thickness_1mm_minSigma.nii.gz").raw
    mask = Bool.(niread("data/AnitaCoverageBrain5/Thickness_1mm_mask.nii.gz").raw)
    χ = niread("data/AnitaCoverageBrain5/Thickness_1mm_MEDI.nii.gz").raw
    
    params = matread("data/AnitaCoverageBrain5/Parameters.mat")

    B0 = params["B0"]
    vsz = tuple(params["Res"]...)
    TEs = params["TE"]
    sz = size(ϕ)
    if !(isdir("./results/$(dataset)"))
        mkdir("./results/$(dataset)")
    end
elseif dataset == "AnitaBrain5b"
    mag = niread("data/AnitaCoverageBrain5/Thickness_1mm_magnitude.nii.gz")
    ϕ = niread("data/AnitaCoverageBrain5/Thickness_1mm_minSigma_vsharp.nii.gz").raw
    mask = Bool.(niread("data/AnitaCoverageBrain5/Thickness_1mm_mask_vsharp.nii.gz").raw)
    χ = niread("data/AnitaCoverageBrain5/Thickness_1mm_MEDI.nii.gz").raw
    
    params = matread("data/AnitaCoverageBrain5/Parameters.mat")

    B0 = params["B0"]
    vsz = tuple(params["Res"]...)
    TEs = params["TE"]
    sz = size(ϕ)
    if !(isdir("./results/$(dataset)"))
        mkdir("./results/$(dataset)")
    end
elseif dataset == "AnitaBrain6"
    mag = niread("data/AnitaCoverageBrain6/Thickness_1mm_magnitude.nii.gz")
    ϕ = niread("data/AnitaCoverageBrain6/Thickness_1mm_minSigma.nii.gz").raw
    mask = Bool.(niread("data/AnitaCoverageBrain6/Thickness_1mm_mask.nii.gz").raw)
    χ = niread("data/AnitaCoverageBrain6/Thickness_1mm_MEDI.nii.gz").raw
    
    params = matread("data/AnitaCoverageBrain5/Parameters.mat")

    B0 = params["B0"]
    vsz = tuple(params["Res"]...)
    TEs = params["TE"]
    sz = size(ϕ)
    if !(isdir("./results/$(dataset)"))
        mkdir("./results/$(dataset)")
    end
elseif dataset == "AnitaBrain6b"
    mag = niread("data/AnitaCoverageBrain6/Thickness_1mm_magnitude.nii.gz")
    ϕ = niread("data/AnitaCoverageBrain6/Thickness_1mm_minSigma_vsharp.nii.gz").raw
    mask = Bool.(niread("data/AnitaCoverageBrain6/Thickness_1mm_mask_vsharp.nii.gz").raw)
    χ = niread("data/AnitaCoverageBrain6/Thickness_1mm_MEDI.nii.gz").raw
    
    params = matread("data/AnitaCoverageBrain5/Parameters.mat")

    B0 = params["B0"]
    vsz = tuple(params["Res"]...)
    TEs = params["TE"]
    sz = size(ϕ)
    if !(isdir("./results/$(dataset)"))
        mkdir("./results/$(dataset)")
    end
elseif dataset == "AnitaBrain7"
    mag = niread("data/AnitaCoverageBrain7/Thickness_1mm_magnitude.nii.gz")
    ϕ = niread("data/AnitaCoverageBrain7/Thickness_1mm_minSigma.nii.gz").raw
    mask = Bool.(niread("data/AnitaCoverageBrain7/Thickness_1mm_mask.nii.gz").raw)
    χ = niread("data/AnitaCoverageBrain7/Thickness_1mm_MEDI.nii.gz").raw
    
    params = matread("data/AnitaCoverageBrain5/Parameters.mat")

    B0 = params["B0"]
    vsz = tuple(params["Res"]...)
    TEs = params["TE"]
    sz = size(ϕ)
    if !(isdir("./results/$(dataset)"))
        mkdir("./results/$(dataset)")
    end
elseif dataset == "AnitaBrain7b"
    mag = niread("data/AnitaCoverageBrain7/Thickness_1mm_magnitude.nii.gz")
    ϕ = niread("data/AnitaCoverageBrain7/Thickness_1mm_minSigma_vsharp.nii.gz").raw
    mask = Bool.(niread("data/AnitaCoverageBrain7/Thickness_1mm_mask_vsharp.nii.gz").raw)
    χ = niread("data/AnitaCoverageBrain7/Thickness_1mm_MEDI.nii.gz").raw
    
    params = matread("data/AnitaCoverageBrain5/Parameters.mat")

    B0 = params["B0"]
    vsz = tuple(params["Res"]...)
    TEs = params["TE"]
    sz = size(ϕ)
    if !(isdir("./results/$(dataset)"))
        mkdir("./results/$(dataset)")
    end
elseif dataset == "AnitaBrain8"
    mag = niread("data/AnitaCoverageBrain8/Thickness_1mm_magnitude.nii.gz")
    ϕ = niread("data/AnitaCoverageBrain8/Thickness_1mm_minSigma.nii.gz").raw
    mask = Bool.(niread("data/AnitaCoverageBrain8/Thickness_1mm_mask.nii.gz").raw)
    χ = niread("data/AnitaCoverageBrain8/Thickness_1mm_MEDI.nii.gz").raw
    
    params = matread("data/AnitaCoverageBrain5/Parameters.mat")

    B0 = params["B0"]
    vsz = tuple(params["Res"]...)
    TEs = params["TE"]
    sz = size(ϕ)
    if !(isdir("./results/$(dataset)"))
        mkdir("./results/$(dataset)")
    end
elseif dataset == "AnitaBrain8b"
    mag = niread("data/AnitaCoverageBrain8/Thickness_1mm_magnitude.nii.gz")
    ϕ = niread("data/AnitaCoverageBrain8/Thickness_1mm_minSigma_vsharp.nii.gz").raw
    mask = Bool.(niread("data/AnitaCoverageBrain8/Thickness_1mm_mask_vsharp.nii.gz").raw)
    χ = niread("data/AnitaCoverageBrain8/Thickness_1mm_MEDI.nii.gz").raw
    
    params = matread("data/AnitaCoverageBrain5/Parameters.mat")

    B0 = params["B0"]
    vsz = tuple(params["Res"]...)
    TEs = params["TE"]
    sz = size(ϕ)
    if !(isdir("./results/$(dataset)"))
        mkdir("./results/$(dataset)")
    end
elseif dataset == "AnitaBrain9"
    mag = niread("data/AnitaCoverageBrain9/Thickness_1mm_magnitude.nii.gz")
    ϕ = niread("data/AnitaCoverageBrain9/Thickness_1mm_minSigma.nii.gz").raw
    mask = Bool.(niread("data/AnitaCoverageBrain9/Thickness_1mm_mask.nii.gz").raw)
    χ = niread("data/AnitaCoverageBrain9/Thickness_1mm_MEDI.nii.gz").raw
    
    params = matread("data/AnitaCoverageBrain5/Parameters.mat")

    B0 = params["B0"]
    vsz = tuple(params["Res"]...)
    TEs = params["TE"]
    sz = size(ϕ)
    if !(isdir("./results/$(dataset)"))
        mkdir("./results/$(dataset)")
    end
elseif dataset == "AnitaBrain9b"
    mag = niread("data/AnitaCoverageBrain9/Thickness_1mm_magnitude.nii.gz")
    ϕ = niread("data/AnitaCoverageBrain9/Thickness_1mm_minSigma_vsharp.nii.gz").raw
    mask = Bool.(niread("data/AnitaCoverageBrain9/Thickness_1mm_mask_vsharp.nii.gz").raw)
    χ = niread("data/AnitaCoverageBrain9/Thickness_1mm_MEDI.nii.gz").raw
    
    params = matread("data/AnitaCoverageBrain5/Parameters.mat")

    B0 = params["B0"]
    vsz = tuple(params["Res"]...)
    TEs = params["TE"]
    sz = size(ϕ)
    if !(isdir("./results/$(dataset)"))
        mkdir("./results/$(dataset)")
    end
end

## ------------------------- ##
## Reconstruction Comparison ##
## ------------------------- ##
algorithms = [ikd, rts, tikh, tsvd, tkd, ndi]
rmse = zeros(length(algorithms))
nrmse = similar(rmse)
xsim = similar(rmse)
psnr = similar(rmse)
nχ = norm(χ)
for (ind, alg) in enumerate(algorithms)
    local savename = "./results/$(dataset)/$(alg)_defaults.nii.gz"
    if isfile(savename)
        @info "Loading $savename"
        rec = niread(savename).raw
    else
        @info "Generating $savename"
        @time rec = alg(ϕ, mask, vsz);
        niwrite(savename,  NIVolume(mag.header, Float32.(rec)))
    end
    rmse[ind] = sqrt(mse(rec, χ))
    nrmse[ind] = rmse[ind]/nχ
    xsim[ind] = assess(XSIM(), mask.*rec, χ)
    psnr[ind] = assess(PSNR(), mask.*rec, χ)
end

# TV "tuned"
Weight = mask.*(mag.raw[:,:,:,1]./maximum(mag.raw))
λTV = 1e-5;
savename = "./results/$(dataset)/nltv_W_lambda_$(λTV).nii.gz"
if isfile(savename)
    recNlTV = niread(savename).raw
else
for retries in 1:10
    global recNlTV = nltv(ϕ, mask, vsz; W=Weight, lambda=λTV, verbose=false)
    if !any(isnan.(recNlTV))
        break
    end
end
niwrite(savename, NIVolume(mag.header, Float32.(recNlTV)))
end
algorithms = cat(algorithms, "nlTV λ$(λTV)"; dims=1)
rmse = cat(rmse, sqrt(mse(mask.*recNlTV, mask.*χ)); dims=1)
nrmse = cat(nrmse, rmse[end]/nχ; dims=1)
xsim = cat(xsim, assess(XSIM(), mask.*recNlTV, mask.*χ); dims=1)
psnr = cat(psnr, assess(PSNR(), mask.*recNlTV, mask.*χ); dims=1)

# CS
D = dipole_kernel(sz, vsz; transform=:fft)
_sz = size(D)
P = plan_fft(ϕ)
iP = inv(P)
δIKD = 0.25
λCS = 1e-5
Sk = abs.(D) .> δIKD
b̂ = ComplexF32.(Sk .* inv.(D) .* (P * (mask .* ϕ)));
niwrite("./results/$(dataset)/b.nii.gz", NIVolume(mag.header, ComplexF32.(b̂)))
niwrite("./results/$(dataset)/Sk.nii.gz", NIVolume(mag.header, Float32.(Sk)))
niwrite("./results/$(dataset)/Schi.nii.gz", NIVolume(mag.header, Float32.(mask)))

function CSfwd!(res, v, α, β)
    res .= vec( (Sk .* (P * ( reshape( v, sz)))))
end
function CSbwd!(res, v, α, β)
    res .= vec( ( (iP * (Sk .* reshape(v, _sz)))))
end
CSop = LinearOperator(eltype(b̂), prod(_sz), prod(sz), false, false,
                    CSfwd!,
                    nothing,
                    CSbwd!);
W = linearOperator("Wavelet", sz); # db2 wavelet
regL1W = Regularization("L1", Float32(λCS); shape=sz, sparseTrafo=W)
solverCS = createLinearSolver("fista", CSop, vec(b̂); reg=regL1W, ρ=1, iterations=50, verbose=false) # has trouble with rfft :-(
savename = "./results/$(dataset)/cs_L1_$(λCS)_delta_$(δIKD).nii.gz"
if isfile(savename)
    recCS = niread(savename).raw
else
recCS = solve(solverCS, vec(b̂))
recCS = real(reshape(recCS,sz))
niwrite(savename, NIVolume(mag.header, Float32.(recCS)))
end
algorithms = cat(algorithms, "cs λ$(λCS) δ$(δIKD)"; dims=1)
rmse = cat(rmse, sqrt(mse(mask.*recCS, mask.*χ)); dims=1)
nrmse = cat(nrmse, rmse[end]/nχ; dims=1)
xsim = cat(xsim, assess(XSIM(), mask.*recCS, mask.*χ); dims=1)
psnr = cat(psnr, assess(PSNR(), mask.*recCS, mask.*χ); dims=1)

# Reg. IKD
function CSIKDfwd!(res, v, α, β)
    res .= vec( (Sk .* (P * (mask .* reshape( v, sz)))))
end
function CSIKDbwd!(res, v, α, β)
    res .= vec( (mask .* (iP * (Sk .* reshape(v, _sz)))))
end
CSIKDop = LinearOperator(eltype(b̂), prod(_sz), prod(sz), false, false,
                    CSIKDfwd!,
                    nothing,
                    CSIKDbwd!);
solverCSIKD = createLinearSolver("fista", CSIKDop, vec(b̂); reg=regL1W, ρ=1, iterations=50, verbose=false) # has trouble with rfft :-(
savename = "./results/$(dataset)/ikdrReg_L1_$(λCS)_delta_$(δIKD).nii.gz"
if isfile(savename)
    recCSIKD = niread(savename).raw
else
    b̂ = ComplexF32.(Sk .* inv.(D) .* (P * (mask .* ϕ)))
    recCSIKD = solve(solverCSIKD, vec(b̂))
    recCSIKD = real(reshape(recCSIKD, sz))
    niwrite(savename, NIVolume(mag.header, Float32.(recCSIKD)))
end
algorithms = cat(algorithms, "reg ikd λ$(λCS) δ$(δIKD)"; dims=1)
rmse = cat(rmse, sqrt(mse(mask.*recCSIKD, mask.*χ)); dims=1)
nrmse = cat(nrmse, rmse[end]/nχ; dims=1)
xsim = cat(xsim, assess(XSIM(), mask.*recCSIKD, mask.*χ); dims=1)
psnr = cat(psnr, assess(PSNR(), mask.*recCSIKD, mask.*χ); dims=1)

# Save results
open("./results/$(dataset)/recon_metrics.csv", "w") do io
    writedlm(io, zip(["Alg";algorithms],["XSIM";xsim],["RMSE";rmse],["NRMSE";nrmse],["PSNR";psnr]), ',')
end

if dataset=="Challenge" && sweepBandLimitFlag
    @info "Running δ parameter sweep"
## ----------------------- ##
## Threshold investigation ##
## ----------------------- ##
thresholds = 0.5:-0.01:0.01;

rmse = zeros(length(thresholds))
xsim = similar(rmse)
psnr = similar(rmse)
for (ind, δ) in enumerate(thresholds)
    local savename = "./results/$(dataset)/ikd_delta_$(δ).nii.gz"
    if isfile(savename)
        @info "Loading $savename"
        recλ = niread(savename).raw
    else
        @info "Generating $savename"
        recλ = ikd(ϕ, mask, vsz; delta=δ)
        niwrite(savename,  NIVolume(mag.header, Float32.(recλ)))
    end
    rmse[ind] = sqrt(mse(recλ, χ))
    xsim[ind] = assess(XSIM(), mask.*recλ, χ)
    psnr[ind] = assess(PSNR(), mask.*recλ, χ)
end
open("./results/$(dataset)/ikd_delta.csv", "w") do io
    writedlm(io, zip(["Delta";thresholds],["XSIM";xsim],["RMSE";rmse],["PSNR";psnr]), ',')
end
end

if dataset=="Challenge" && sweepMaskFlag
    @info "Running Sχ parameter sweep"
## -------------------##
## Mask investigation ##
## -------------------##
dilations = 1:8;
erosions = 1:8;
masksizes = [-collect(erosions)...,collect(dilations)...,0]

Sχ = similar(mask);

rmse = zeros(length(masksizes))
xsim = similar(rmse)
psnr = similar(rmse)
copyto!(Sχ, mask);
for (ind, vxnum) in enumerate(erosions)
    Sχ = erode(Sχ);
    local savename = "./results/$(dataset)/ikd_mask_$(vxnum)_eroded.nii.gz"
    if isfile(savename)
        @info "Loading $savename"
        recErode = niread(savename).raw
    else
        @info "Generating $savename"
        recErode = ikd(ϕ.*Sχ, Sχ, vsz; tol=1e-2);
        niwrite(savename,  NIVolume(mag.header, Float32.(recErode)));
        niwrite("./results/$(dataset)/mask_$(vxnum)_eroded.nii.gz", NIVolume(mag.header, Float32.(Sχ)))
    end
    rmse[ind] = sqrt(mse(recErode, χ))
    xsim[ind] = assess(XSIM(), recErode, χ)
    psnr[ind] = assess(PSNR(), recErode, χ)
end

copyto!(Sχ, mask);
for (ind, vxnum) in enumerate(dilations)
    Sχ = dilate(Sχ);
    local savename = "./results/$(dataset)/ikd_mask_$(vxnum)_dilated.nii.gz"
    if isfile(savename)
        @info "Loading $savename"
        recDilate = niread(savename).raw
    else
        @info "Generating $savename"
        recDilate = ikd(ϕ.*Sχ, Sχ, vsz; tol=1e-2);
        niwrite(savename,  NIVolume(mag.header, Float32.(recDilate)));
        niwrite("./results/$(dataset)/mask_$(vxnum)_dilated.nii.gz", NIVolume(mag.header, Float32.(Sχ)))
    end
    rmse[ind+length(erosions)] = sqrt(mse(recDilate, χ))
    xsim[ind+length(erosions)] = assess(XSIM(), recDilate, χ)
    psnr[ind+length(erosions)] = assess(PSNR(), recDilate, χ)
end

copyto!(Sχ, mask);
savename = "./results/$(dataset)/ikd_mask_0.nii.gz"
if isfile(savename)
    @info "Loading $savename"
    recOrigMask = niread(savename).raw
else
    @info "Generating $savename"
    recOrigMask = ikd(ϕ.*Sχ, Sχ, vsz; tol=0.1)
    niwrite(savename,  NIVolume(mag.header, Float32.(recOrigMask)))
    niwrite("./results/$(dataset)/mask_0.nii.gz", NIVolume(mag.header, Float32.(Sχ)))
end
rmse[end] = sqrt(mse(recOrigMask, χ))
xsim[end] = assess(XSIM(), recOrigMask, χ)
psnr[end] = assess(PSNR(), recOrigMask, χ)

open("./results/$(dataset)/ikd_mask.csv", "w") do io
    writedlm(io, zip(["Mask";masksizes],["XSIM";xsim],["RMSE";rmse],["PSNR";psnr]), ',')
end

end


if dataset=="Challenge" && tuningFlag
    ## ----------------------- ##
    ## Parameter investigation ##
    ## ----------------------- ##
    # FANSI
    @info "Running λTV parameter sweep"
    thresholds = vcat([10f1 .^ (0:-0.5:-4), collect(0.008:-0.002:0.0001), collect(0.0008:-0.0002:0.00001)]...)
    Weight = mask.*(mag.raw[:,:,:,1]./maximum(mag.raw))

    rmse = zeros(length(thresholds))
    xsim = similar(rmse)
    psnr = similar(rmse)
    for (ind, λTV) in enumerate(thresholds)
        local savename = "./results/$(dataset)/nlTV_lambda_$(λTV).nii.gz"
        if isfile(savename)
            recλ = niread(savename).raw
        else
            for retries in 1:10
                recλ = nltv(ϕ, mask, vsz; W=Weight, lambda=λTV, verbose=false)
                if !any(isnan.(recλ))
                    niwrite(savename,  NIVolume(mag.header, Float32.(recλ)));
                    break
                end
                
            end
        end
        rmse[ind] = sqrt(mse(mask.*recλ, mask.*χ))
        xsim[ind] = assess(XSIM(), mask.*recλ, mask.*χ)
        psnr[ind] = assess(PSNR(), mask.*recλ, mask.*χ)
    end
    open("./results/$(dataset)/nlTV_lambda.csv", "w") do io
        writedlm(io, zip(["LambdaTV";thresholds],["XSIM";xsim],["RMSE";rmse],["PSNR";psnr]), ',')
    end

    # CS
    @info "Running λCS parameter sweep"
    D = dipole_kernel(sz, vsz; transform=:fft)
    _sz = size(D)
    P = plan_fft(ϕ)
    iP = inv(P)
    δIKD = 0.28
    λCS = 0.001
    Sk = abs.(D) .> δIKD

    function CSfwd!(res, v, α, β)
        res .= vec( (Sk .* (P * ( reshape( v, sz)))))
    end
    function CSbwd!(res, v, α, β)
        res .= vec( ( (iP * (Sk .* reshape(v, _sz)))))
    end
    CSop = LinearOperator(eltype(b̂), prod(_sz), prod(sz), false, false,
                        CSfwd!,
                        nothing,
                        CSbwd!)
    W = linearOperator("Wavelet", sz)

    thresholds = vcat(sort([1f2 .^ (0:-0.5:-2); 0.5*(1f2 .^ (0:-0.5:-2))]), 1e-5, 5e-6, 1e-6, 5e-7)

    rmse = zeros(length(thresholds))
    xsim = similar(rmse)
    psnr = similar(rmse)
    for (ind, λCS) in enumerate(thresholds)
        local savename = "./results/$(dataset)/cs_lambda_$(λCS).nii.gz"
        if isfile(savename)
            recλ = niread(savename).raw
        else
            local b̂ = ComplexF32.(Sk .* inv.(D) .* (P * (mask .* ϕ)))
            local regL1W = Regularization("L1", Float32(λCS); shape=sz, sparseTrafo=W)
            local solverCS = createLinearSolver("fista", CSop, vec(b̂); reg=regL1W, ρ=1, iterations=50, verbose=false)
            
            local b̂ = ComplexF32.(Sk .* inv.(D) .* (P * (mask .* ϕ)))
            recλ = solve(solverCS, vec(b̂))
            recλ = real(reshape(recλ,sz))
            niwrite(savename,  NIVolume(mag.header, Float32.(recλ)))
        end
        rmse[ind] = sqrt(mse(mask.*recλ, mask.*χ))
        xsim[ind] = assess(XSIM(), mask.*recλ, mask.*χ)
        psnr[ind] = assess(PSNR(), mask.*recλ, mask.*χ)
    end
    open("./results/$(dataset)/cs_lambda.csv", "w") do io
        writedlm(io, zip(["LambdaCS";thresholds],["XSIM";xsim],["RMSE";rmse],["PSNR";psnr]), ',')
    end

end
