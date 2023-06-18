# cameramanFT=fft(Float64.(cameraman))
# dct_fun = plan_dct(cameraman,dims=512)
# @view(abs.(cameramanFT[1:end,1:end]))
# barbara = testimage("barbara_gray_512.bmp")
# lena = testimage("lena_gray_512.tif")
# chart = testimage("resolution_test_512.tif")
# x = vec(cameraman);
# rng = MersenneTwister(20230511)
# # n = length(x)
# # m = 100
# # r = 64
# # A = rand(rng, m, n)
# # w = rand(rng, m)
# # b = A * x + w |> y->Float64.(y)
# # ATA = A'*A
# m, n, r = 512, 2048, 64
# x = zeros(n)
# q = randperm(n)
# A = randn(rng, m, n)
# B = qr(A).Q * Matrix(I, size(A)...)
# b = A * x + vec(0.001 * rand(rng, m, 1)) #|> y->Float64.(x)
# plot(b)
# denoise(b) |> plot

# # img = load("./imgs/lenna.png")
# # ndims(img)
# # x = permutedims(Float64.(lena))
# # L = 2
# # xts = wplotim(x, L, wavelet(WT.db3, WT.Filter))
# # Gray.(xts)

include("dependencies.jl")
using Images, TestImages, FileIO
using Wavelets
using FFTW
# FFTW.set_provider!("mkl")
# using ProgressBars

source_image = load("./imgs/man.bmp");
source_image_matrix = Float64.(Gray.(source_image));
(m, n) = size(source_image_matrix)
middle = fld(n, 2) + 1

sigma = sqrt(0.01)
h = zeros(m, n)
for i = -4:4
    for j = -4:4
        h[i+middle, j+middle] = (1 / (1 + i * i + j * j))
    end
end
h = fftshift(h)
h = h / sum(h)
R(x) = real(ifft(fft(h) .* fft(x)))
# real(ifft2(conj(fft2(h)).*fft2(x)));
RT(x) = real(ifft(conj(fft(h)) .* fft(x)))
mywt = wavelet(WT.db2)
W(x) = idwt(x, mywt, 3)
myWT(x) = dwt(x, mywt, 3)
# W = @(x) midwt(x,wav,3);
# WT = @(x) mdwt(x,wav,3);
y = R(source_image_matrix) + sigma * randn(m, n);
mosaicview(Gray.(y), source_image; nrow=1, npad=20, fillvalue=0.5)

Afun(x) = R(W(x))
ATfun(x) = myWT(RT(x))
x0 = ATfun(y)
τ = 0.35

Aty = ATfun(y)

# % Initialization
x = x0
true_x = myWT(source_image_matrix)
# Gray.(true_x)
# max_tau = norm(x0, Inf)
u = max.(x, 0)
v = max.(-x, 0)
num_nz_x = count(x -> x != 0, x)
final_tau = τ
final_stopCriterion = 1
final_tolA = 1.e-5
cont_factors = 1;
cont_steps = 1;
mses_matrix = Vector{Float64}(undef, 2000)
iter = 1
mses_matrix[iter] = norm(x - true_x)^2
resid = y - Afun(x)

c = τ * ones(2n) + vcat(-x0, x0)
z0 = vcat(max.(x0, 0), max.(-x0, 0))
original_image_signal = myWT(source_image_matrix)

