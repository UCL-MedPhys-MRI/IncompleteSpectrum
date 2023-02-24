using DelimitedFiles
using LaTeXStrings
using MIRTjim: jim, mid3
using NIfTI: niread
using Plots
using QSM: padarray!

gr()

function plotslices(x;kwargs...)
    x = Array(x)
    sz = size(x)
    mid_slice = sz .÷ 2
    padsize = (maximum(sz), maximum(sz))
    if eltype(x) == Bool
        y = ones(sz)
        y[.!x] .= 0
        x = y
    end
    axial = Array{eltype(x)}(undef,padsize)
    coronal = similar(axial)
    sagittal = similar(coronal)
    fillel = NaN    
    padarray!(axial, x[:,:,mid_slice[3]], :fill, fillel)
    padarray!(coronal, x[:,mid_slice[2],:], :fill, fillel)
    padarray!(sagittal, x[mid_slice[1],:,:], :fill, fillel)
    planes = vcat(axial, coronal, sagittal)
    jim(planes; size = (900, 350),
                axis=true, 
                colorbar=true, 
                yflip=false,
                background_color = :transparent,
                xticks=[], yticks=[],
                xlims=(1,size(planes)[1]), ylims=(1,size(planes)[2]),
                kwargs...)
end

x = niread("./data/ChallengeSim2SNR1/Chi.nii.gz");
support = niread("./data/ChallengeSim2SNR1/Brain.nii.gz");

Plots.gr_cbar_width[] = 0.02
pl1 = plotslices(x.raw./support.raw;
        clim=(-0.1,0.1),
        colorbar_title="ppm", 
        right_margin = 20Plots.px, 
        fg = "black",
        size = (900, 300), )
savefig("./figures/groundtruth.png")

function plotaxial(x; kwargs...)

    sz = size(x)
    mid_slice = sz .÷ 2

    axial = x[:,:,mid_slice[3],:]
    axial = reshape(axial, (sz[1], :))
    jim(axial;
    size = (1200, 800), 
    axis=true, 
    colorbar=true, 
    yflip=false,
    xticks=[], yticks=[],
    xlims=(1,size(axial)[1]), ylims=(1,size(axial)[2]),
    kwargs...)
end

function plotsagittal(x; kwargs...)

    sz = size(x)
    mid_slice = sz .÷ 2

    sagittal = x[mid_slice[1],:,:,:]
    sagittal = reshape(sagittal, (sz[2], :))
    jim(sagittal;
    size = (1200, 800), 
    axis=true, 
    colorbar=true, 
    yflip=false,
    xticks=[], yticks=[],
    xlims=(1,size(sagittal)[1]), ylims=(1,size(sagittal)[2]),
    kwargs...)
end

stats = readdlm("./results/Challenge/ikd_delta.csv", ',')

λs = stats[18:end-1,1]
psnrs = stats[18:end-1,4]
xsims = stats[18:end-1,2]
optλ = [findmax(x) for x in (psnrs,xsims)]

cur_colors = theme_palette(:auto);
plt1a = plot(λs, psnrs; 
        lw = 3, 
        c = cur_colors[1],
        ylabel = "PSNR", 
        xlabel = L"t_{well}"*"\n(a)",
        legend = false);
scatter!((λs[optλ[1][2]],optλ[1][1]); marker = (:diamond, 5, 5, cur_colors[4]))

plt1b = plot(λs, xsims;
        lw = 3, 
        c = cur_colors[2],
        ylabel = "XSIM", 
        xlabel = L"t_{well}"*"\n(c)",
        legend = false);
scatter!((λs[optλ[2][2]],optλ[2][1]); marker = (:diamond, 5, 5, cur_colors[4]))

χikd_PSRN = niread("./results/Challenge/ikd_delta_$(λs[optλ[1][2]]).nii.gz");
χikd_XSIM = niread("./results/Challenge/ikd_delta_$(λs[optλ[2][2]]).nii.gz");
support = niread("./data/ChallengeSim2SNR1/Brain.nii.gz");
bbox = ((20,14,  21),
        (147,183,170));
χikd_PSRN = χikd_PSRN.raw[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]];
χikd_XSIM = χikd_XSIM.raw[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]];
support = support.raw[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]];

plt2a = plotsagittal(χikd_PSRN./support;clim=(-0.1,0.1),xlabel="\n(b)",colorbar_title="ppm");
quiver!([135+17], [72-15], quiver=([-17], [15]);arrow=arrow(), linewidth=2, c=:orange);
plt2b = plotsagittal(χikd_XSIM./support;clim=(-0.1,0.1),xlabel="\n(d)",colorbar_title="ppm");
quiver!([135+17], [72-15], quiver=([-17], [15]);arrow=arrow(), linewidth=2, c=:orange);

jim(plt1a,plt2a,plt1b,plt2b; 
                left_margin = 18Plots.px,
                size = (900, 600), 
                background_color="transparent")
savefig("./figures/psnr_v_xsim_full.png")

GT = niread("./data/ChallengeSim2SNR1/Chi.nii.gz");
GT = GT.raw[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]];

plt5a = plotsagittal((χikd_PSRN .- χikd_XSIM)./support;clim=(-0.1,0.1),xlabel="\n(a)",colorbar_title="ppm");
plt5b = plotsagittal((χikd_PSRN .- GT)./support;clim=(-0.1,0.1),xlabel="\n(b)",colorbar_title="ppm");
plt5c = plotsagittal((χikd_XSIM .- GT)./support;clim=(-0.1,0.1),xlabel="\n(c)",colorbar_title="ppm");
jim(plt5a,plt5b,plt5c; 
                bottom_margin = 18Plots.px,
                size = (600, 1200), 
                background_color="transparent",layout=(3,1))
savefig("./figures/psnr_v_xsim_diffs_gt.png")
jim(plt5a;
                size = (600, 400), 
                background_color="transparent",layout=(3,1))
savefig("./figures/psnr_min_xsim.png")

## Figure 2
stats = readdlm("./results/Challenge/ikd_mask.csv", ',')
masksizes = stats[2:end,1]
psnrs = stats[2:end,4]
xsims = stats[2:end,2]

I = sortperm(masksizes)
plt1 = plot(masksizes[I], psnrs[I];
        lw=2, c=cur_colors[1],
        ylabel="PSNR", xlabel="(a)",
        marker = :utriangle,
        markersize = 6,
        markerstrokecolor = cur_colors[1],
        legend = false);
plt1b = plot(masksizes[I], xsims[I];
        lw=2, c=cur_colors[1],
        ylabel="XSIM", xlabel="Number of voxels of mask erosion / dilation\n(b)",
        marker = :utriangle,
        markersize = 6,
        markerstrokecolor = cur_colors[1],
        legend = false);

χikd_dilated = niread("./results/Challenge/ikd_mask_5_dilated.nii.gz");
support_dilated = niread("./results/Challenge/mask_5_dilated.nii.gz");
χikd_eroded = niread("./results/Challenge/ikd_mask_5_eroded.nii.gz");
support_eroded  = niread("./results/Challenge/mask_5_eroded.nii.gz");
bbox = ((15,9,  19),
        (153,188,175));
χikd_dilated = χikd_dilated.raw[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]];
χikd_eroded = χikd_eroded.raw[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]];
support_dilated = support_dilated.raw[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]];
support_eroded = support_eroded.raw[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]];

plt2a = plotsagittal(χikd_eroded./support_eroded;clim=(-0.1,0.1),xlabel="(c)",colorbar_title="ppm");
plt2b = plotsagittal(χikd_dilated./support_dilated;clim=(-0.1,0.1),xlabel="(d)",colorbar_title="ppm");

jim(jim(plt1,plt1b;layout=(2,1)), jim(plt2a, plt2b;layout=(2,1));
        left_margin = 18Plots.px,
        bottom_margin = 10Plots.px, 
        size = (900, 600), 
        background_color="white")

 savefig("./figures/mask_erosion_dilation_psnrxsim.png")

#### Figure 3
using QSM: dipole_kernel

recon = niread("./results/Challenge/ikd_defaults.nii.gz").raw;
band_limit = niread("./results/Challenge/Sk.nii.gz");
x = niread("./data/ChallengeSim2SNR1/Chi.nii.gz");
b = niread("./results/Challenge/b.nii.gz").raw;
b[isinf.(b)] .= 0;

sz = size(b)
D = dipole_kernel(sz, (1,1,1); method=:k)
Sk  = abs.(D) .> 0.25
using FFTW: plan_fft, fftshift
F = plan_fft(x)
recon = F*recon;

plt = plotsagittal(cat(fftshift(recon),fftshift(b);dims=4);clim=(0,100))

plt1 = plotsagittal(fftshift(b);clim=(0,100),colorbar=false,xlabel="(a)")
plt2 = plotsagittal(fftshift(Sk);clim=(0,1),colorbar=false,xlabel="(b)")
annotate!(90,180, ("Well posed", :k, :center, 8))
annotate!(90,168, (L"(|D|>t_{well})", :k, :center, 8))
annotate!(30,160, ("Ill posed", :white, :center, 8))
annotate!(30,148, (L"(|D|<t_{ill})", :white, :center, 8))

plt3 = plotsagittal(fftshift(recon);clim=(0,100),colorbar=false, xlabel="(c)")
jim(plt1,plt2,plt3;layout=(1,3),size = (900, 300), bottom_margin = 15Plots.px)

savefig("./figures/splitspectrum.png")

### Figure 5
recCS = niread("./results/AnitaBrain9/cs_L1_5.0e-6_delta_0.28.nii.gz").raw;
recIKD = niread("./results/AnitaBrain9/ikd_defaults.nii.gz").raw;
recIKDCS = niread("./results/AnitaBrain9/ikdrReg_L1_5.0e-6_delta_0.28.nii.gz").raw;
recNlTV = niread("./results/AnitaBrain9/nltv_W_lambda_1.0e-5.nii.gz").raw;
recTKD = niread("./results/AnitaBrain9/tkd_defaults.nii.gz").raw;
recNDI = niread("./results/AnitaBrain9/ndi_defaults.nii.gz").raw;
support = niread("./data/AnitaCoverageBrain9/Thickness_1mm_mask.nii.gz").raw;

bbox = ((49,32,  1),
        (186,203,140));
recIKD = recIKD[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]];
support = support[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]];
recCS = recCS[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]]
recIKDCS = recIKDCS[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]]
recNlTV = recNlTV[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]]
recTKD = recTKD[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]]
recNDI = recNDI[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]]

function plotsagittal(x; kwargs...)

    sz = size(x)
    mid_slice = sz .÷ 2

    sagittal = x[mid_slice[1]+15,:,:,:]
    sagittal = reshape(sagittal, (sz[2], :))
    jim(sagittal;
    size = (1200, 800), 
    axis=true, 
    colorbar=true, 
    yflip=false,
    background_color = :transparent,
    xticks=[], yticks=[],
    xlims=(1,size(sagittal)[1]), ylims=(1,size(sagittal)[2]),
    kwargs...)
end

function plotcoronal(x; kwargs...)

        sz = size(x)
        mid_slice = sz .÷ 2
    
        coronal = x[:,mid_slice[2]+4,24:end,:]
        coronal = reshape(coronal, (sz[1], :))
        jim(coronal;
        size = (1200, 800), 
        axis=true, 
        colorbar=true, 
        yflip=false,
        background_color = :transparent,
        xticks=[], yticks=[],
        xlims=(1,size(coronal)[1]), ylims=(1,size(coronal)[2]),
        kwargs...)
end

Plots.gr_cbar_width[] = 0.02
l = @layout [a{0.95w} b]
cb = heatmap([NaN;;], framestyle=:none, clims=(-0.1,0.1), cbar=true, cmap=:grays,size=(50,300),colorbar_title="ppm")

plt1a = plotsagittal(recIKD./support     ; xlabel="Incomplete Spectrum\n(a)", clim=(-0.1,0.1));
quiver!([39-17], [22-15], quiver=([17], [15]);arrow=arrow(), linewidth=2, c=:orange);
plt2a = plotsagittal(recIKDCS./support   ; xlabel="IS Regularised\n(b)", clim=(-0.1,0.1));
quiver!([39-17], [22-15], quiver=([17], [15]);arrow=arrow(), linewidth=2, c=:orange);
plt3a = plotsagittal(recCS./support      ; xlabel="Compressed Sensing\n(c)", clim=(-0.1,0.1));
plt4a = plotsagittal(recTKD./support     ; xlabel="TKD\n(d)", clim=(-0.1,0.1));
plt5a = plotsagittal(recNlTV./support    ; xlabel="FANSI\n(e)", clim=(-0.1,0.1));
plt6a = plotsagittal(recNDI./support     ; xlabel="NDI\n(f)", clim=(-0.1,0.1));

plt1 = jim(plt1a, plt2a, plt3a, plt4a, plt5a, plt6a;
        fg="black", size = (900, 500), clim=(-0.1,0.1), colorbar=false,
        bottom_margin = 25Plots.px,)
plt1 = jim(plt1,cb, layout=l)

plt1b = plotcoronal(recIKD./support     ; xlabel="Incomplete Spectrum\n(g)", clim=(-0.1,0.1))
plt2b = plotcoronal(recIKDCS./support   ; xlabel="IS Regularised\n(h)", clim=(-0.1,0.1));
plt3b = plotcoronal(recCS./support      ; xlabel="Compressed Sensing\n(i)", clim=(-0.1,0.1));
plt4b = plotcoronal(recTKD./support     ; xlabel="TKD\n(j)", clim=(-0.1,0.1));
plt5b = plotcoronal(recNlTV./support    ; xlabel="FANSI\n(k)", clim=(-0.1,0.1));
plt6b = plotcoronal(recNDI./support     ; xlabel="NDI\n(l)", clim=(-0.1,0.1));

plt2 = jim(plt1b, plt2b, plt3b, plt4b, plt5b, plt6b;
        fg="black", size = (900, 500), clim=(-0.1,0.1), colorbar=false,
        bottom_margin = 25Plots.px,)
plt2 = jim(plt2,cb, layout=l)

jim(plt1,plt2;layout=(2,1), 
        size = (900, 1100), 
        background_color="transparent")

savefig("./figures/anita_comparison.png")

## Differences

plt1a = plotsagittal((recIKD  .- recIKDCS)./support   ; xlabel="IS - IS Reg\n(a)", clim=(-0.1,0.1))
plt2a = plotsagittal((recIKDCS .- recCS)./support    ; xlabel="IS Reg - CS\n(b)", clim=(-0.1,0.1));
plt3a = plotsagittal((recIKD .- recCS)./support      ; xlabel="IS - CS\n(c)", clim=(-0.1,0.1));
plt4a = plotsagittal((recIKD .- recTKD)./support     ; xlabel="IS - TKD\n(d)", clim=(-0.1,0.1));
plt5a = plotsagittal((recIKD .- recNlTV)./support    ; xlabel="IS - FANSI\n(e)", clim=(-0.1,0.1));
plt6a = plotsagittal((recIKD .- recNDI)./support     ; xlabel="IS - NDI\n(f)", clim=(-0.1,0.1));
quiver!([105+20], [65-15], quiver=([-20], [15]);arrow=arrow(), linewidth=2, c=:orange);

plt1 = jim(plt1a, plt2a, plt3a, plt4a, plt5a, plt6a;
        fg="black", size = (900, 500), clim=(-0.1,0.1), colorbar=false,
        bottom_margin = 25Plots.px,)
plt1 = jim(plt1,cb, layout=l)

plt1b = plotcoronal((recIKD  .- recIKDCS)./support   ; xlabel="IS - IS Reg\n(g)", clim=(-0.1,0.1))
plt2b = plotcoronal((recIKDCS .- recCS)./support    ; xlabel="IS Reg - CS\n(h)", clim=(-0.1,0.1));
plt3b = plotcoronal((recIKD .- recCS)./support      ; xlabel="IS - CS\n(i)", clim=(-0.1,0.1));
plt4b = plotcoronal((recIKD .- recTKD)./support     ; xlabel="IS - TKD\n(j)", clim=(-0.1,0.1));
plt5b = plotcoronal((recIKD .- recNlTV)./support    ; xlabel="IS - FANSI\n(k)", clim=(-0.1,0.1));
plt6b = plotcoronal((recIKD .- recNDI)./support     ; xlabel="IS - NDI\n(l)", clim=(-0.1,0.1));
quiver!([50-17], [40-15], quiver=([17], [15]);arrow=arrow(), linewidth=2, c=:orange);
quiver!([95+17], [40-15], quiver=([-17], [15]);arrow=arrow(), linewidth=2, c=:orange);

plt2 = jim(plt1b, plt2b, plt3b, plt4b, plt5b, plt6b;
        fg="black", size = (900, 500), clim=(-0.1,0.1), colorbar=false,
        bottom_margin = 25Plots.px,)
plt2 = jim(plt2,cb, layout=l)

jim(plt1,plt2;layout=(2,1), 
        size = (900, 1100), 
        background_color="transparent")

savefig("./figures/anita_differences.png")


## Supplemental figures
recCS = niread("./results/AnitaBrain9/cs_L1_0.01_delta_0.28.nii.gz").raw;
recIKD = niread("./results/AnitaBrain9/ikd_defaults.nii.gz").raw;
recIKDCS = niread("./results/AnitaBrain9/ikdrReg_L1_0.01_delta_0.28.nii.gz").raw;
support = niread("./data/AnitaCoverageBrain9/Thickness_1mm_mask.nii.gz").raw;

bbox = ((49,32,  1),
        (186,203,140));
recIKD = recIKD[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]];
support = support[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]];
recIKDCS = recIKDCS[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]];
recCS = recCS[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]];
        
Plots.gr_cbar_width[] = 0.02
l = @layout [a{0.95w} b]
cb = heatmap([NaN;;], framestyle=:none, clims=(-0.1,0.1), cbar=true, cmap=:grays,size=(50,300),colorbar_title="ppm")

plta = plotsagittal(recIKD./support     ; xlabel="Incompl Spec\n(c)", clim=(-0.1,0.1));
pltb = plotsagittal(recIKDCS./support     ; xlabel="Incompl Spec Regularised\n(a)", clim=(-0.1,0.1));
pltc = plotsagittal(recCS./support     ; xlabel="Compressed Sensing\n(b)", clim=(-0.1,0.1));
pltd = plotsagittal((recIKDCS .- recCS)./support     ; xlabel="IS Reg - CS\n(d)", clim=(-0.1,0.1));

plt1 = jim(pltb, pltc, plta;
        fg="black", 
        size = (900, 300),
        layout=(1,3),
        clim=(-0.1,0.1), 
        colorbar=false,
        bottom_margin = 25Plots.px,
        )
plt2 = jim(plt1,cb, layout=l)
savefig("./figures/anita_CS_l1_0.01_comparison.png")

pltd = plotsagittal((recIKDCS .- recCS)./support     ; 
xlabel="IS Reg - CS", 
clim=(-0.05,0.05),
fg="black",
size = (400,300),
colorbar_title="ppm",
right_margin = 5Plots.px,
)
savefig("./figures/anita_CS_l1_0.01_difference.png")


## Commenter figures BGFR
recPDF = niread("./results/AnitaBrain6/ikd_defaults.nii.gz").raw;
recVSHARP = niread("./results/AnitaBrain6b/ikd_defaults.nii.gz").raw;
support = niread("./data/AnitaCoverageBrain6/Thickness_1mm_mask.nii.gz").raw;

bbox = ((55,15,  2),
        (188,182,138));
recPDF = recPDF[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]];
support = support[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]];
recVSHARP = recVSHARP[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]];


plta = plotsagittal(recPDF./support     ; xlabel="PDF\n(a)", clim=(-0.15,0.15));
pltb = plotsagittal(recVSHARP./support     ; xlabel="VSHARP\n(b)", clim=(-0.15,0.15));
pltc = plotaxial(recPDF./support     ; xlabel="PDF\n(a)", clim=(-0.15,0.15));
pltd = plotaxial(recVSHARP./support     ; xlabel="VSHARP\n(b)", clim=(-0.15,0.15));

plt1 = jim(plta, pltb, pltc,pltd;
        fg="black", 
        size = (600, 600),
        layout=(2,2),
        )
savefig("./figures/anita_6_vsharp_v_pdf.png")


### Figure 5
recCS = niread("./results/Challenge/cs_L1_5.0e-6_delta_0.28.nii.gz").raw;
recIKD = niread("./results/Challenge/ikd_defaults.nii.gz").raw;
recIKDCS = niread("./results/Challenge/ikdrReg_L1_5.0e-6_delta_0.28.nii.gz").raw;
recNlTV = niread("./results/Challenge/nltv_W_lambda_1.0e-5.nii.gz").raw;
recTKD = niread("./results/Challenge/tkd_defaults.nii.gz").raw;
recNDI = niread("./results/Challenge/ndi_defaults.nii.gz").raw;
support = niread("./data/ChallengeSim2SNR1/Brain.nii.gz").raw;

bbox = ((19,17,  18),
        (150,184,170));
recIKD = recIKD[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]];
support = support[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]];
recCS = recCS[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]]
recIKDCS = recIKDCS[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]]
recNlTV = recNlTV[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]]
recTKD = recTKD[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]]
recNDI = recNDI[bbox[1][1]:bbox[2][1],bbox[1][2]:bbox[2][2],bbox[1][3]:bbox[2][3]]


Plots.gr_cbar_width[] = 0.02
l = @layout [a{0.95w} b]
cb = heatmap([NaN;;], framestyle=:none, clims=(-0.1,0.1), cbar=true, cmap=:grays,size=(50,300),colorbar_title="ppm")

plt1a = plotsagittal(recIKD./support     ; xlabel="Incomplete Spectrum\n(a)", clim=(-0.1,0.1));
plt2a = plotsagittal(recIKDCS./support   ; xlabel="IS Regularised\n(b)", clim=(-0.1,0.1));
plt3a = plotsagittal(recCS./support      ; xlabel="Compressed Sensing\n(c)", clim=(-0.1,0.1));
plt4a = plotsagittal(recTKD./support     ; xlabel="TKD\n(d)", clim=(-0.1,0.1));
plt5a = plotsagittal(recNlTV./support    ; xlabel="FANSI\n(e)", clim=(-0.1,0.1));
plt6a = plotsagittal(recNDI./support     ; xlabel="NDI\n(f)", clim=(-0.1,0.1));

plt1 = jim(plt1a, plt2a, plt3a, plt4a, plt5a, plt6a;
        fg="black", size = (900, 500), clim=(-0.1,0.1), colorbar=false,
        bottom_margin = 25Plots.px,)
plt1 = jim(plt1,cb, layout=l)

plt1b = plotcoronal(recIKD./support     ; xlabel="Incomplete Spectrum\n(g)", clim=(-0.1,0.1))
plt2b = plotcoronal(recIKDCS./support   ; xlabel="IS Regularised\n(h)", clim=(-0.1,0.1));
plt3b = plotcoronal(recCS./support      ; xlabel="Compressed Sensing\n(i)", clim=(-0.1,0.1));
plt4b = plotcoronal(recTKD./support     ; xlabel="TKD\n(j)", clim=(-0.1,0.1));
plt5b = plotcoronal(recNlTV./support    ; xlabel="FANSI\n(k)", clim=(-0.1,0.1));
plt6b = plotcoronal(recNDI./support     ; xlabel="NDI\n(l)", clim=(-0.1,0.1));

plt2 = jim(plt1b, plt2b, plt3b, plt4b, plt5b, plt6b;
        fg="black", size = (900, 500), clim=(-0.1,0.1), colorbar=false,
        bottom_margin = 25Plots.px,)
plt2 = jim(plt2,cb, layout=l)

jim(plt1,plt2;layout=(2,1), 
        size = (900, 1100), 
        background_color="white")

savefig("./figures/challenge_comparison.png")