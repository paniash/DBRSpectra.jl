using LinearAlgebra
using PyPlot

#%%
function refractive_and_depths(N, n1, n2, d1, d2)
    n = Complex{Float64}[]  # To allow for complex refractive indices
    d = Float64[]
    for i in 1:N
        if i % 2 != 0
            push!(n, n1)
            push!(d, d1)
        else
            push!(n, n2)
            push!(d, d2)
        end
    end

    return n, d
end

function find_thetas(θi, n_med, n_list)
    n0_sinθ_0 = n_med * sin(θi)
    θ_list = Complex{Float64}[]
    l = length(n_list)

    for i in 1:l
        n = n_list[i]
        n_real = real(n)
        n_imag = imag(n)
        θ = asin(n0_sinθ_0/n)
        n_cosθ = sqrt(n^2 - n0_sinθ_0^2)

        check_absorption = n_real * n_imag
        if real(n_cosθ) < -100 * eps(Float64)
            θ = π - θ
        end

        if check_absorption <= abs(100 * eps(Float64))
            if n_real > n0_sinθ_0
                if imag(n_cosθ) < -100 * eps(Float64)
                    θ = π - θ
                end

            else
                if real(n_cosθ) < -100 * eps(Float64)
                    θ = π - θ
                end
            end
        end

        push!(θ_list, θ)
    end

    return θ_list
end

function phase_of_reflectivity(r)
    re = real(r)
    im = imag(r)
    return atan(re/im)
end

function r12_s(n1, n2, θ1, θ2)
    num = n1 * cos(θ1) - n2 * cos(θ2)
    denom = n1 * cos(θ1) + n2 * cos(θ2)
    return num/denom
end

function t12_s(n1, n2, θ1, θ2)
    num = 2 * n1 * cos(θ1)
    denom = n1 * cos(θ1) + n2 * cos(θ2)
    return num/denom
end

function r12_p(n1, n2, θ1, θ2)
    num = n2 * cos(θ1) - n1 * cos(θ2)
    denom = n2 * cos(θ1) + n1 * cos(θ2)
    return num/denom
end

function t12_p(n1, n2, θ1, θ2)
    num = 2 * n1 * cos(θ1)
    denom = n2 * cos(θ1) + n1 * cos(θ2)
    return num/denom
end
#%%

function refractive_metal(λ, λp, λc)
    ωp = 1.0 / λp
    ω = 1.0 / λ
    ωc = 1.0 / λc
    num = ωp^2
    denom = ω^2 + (1.0im) * ω * ωc
    nsq = 1 - num/denom
    n = sqrt(nsq)
    return n
end

function transfer_matrix(wavelength, polarization, θi, n_med, n_list, d_list, metal::Vector{Any})
    if metal[1] == true
        λp = metal[2]
        λc = metal[3]
        metal_depth = metal[4]
        n_metal = refractive_metal(wavelength, λp, λc)

        if metal[5] == "right"
            push!(n_list, n_metal, 1)
            push!(d_list, metal_depth, 10)

        elseif metal[5] == "left"
            # n_list = vcat(n_metal, n_list, 1)
            # d_list = vcat(metal_depth, d_list, 10)
            n_list = vcat(n_metal, n_list[2], n_list, 1)
            d_list = vcat(metal_depth, d_list[2] - metal_depth, d_list, 10)
        end

    else
        push!(n_list, 1)
        push!(d_list, 10)
    end

    l = length(n_list)

    θ_list = find_thetas(θi, n_med, n_list)
    r_nnp1 = Complex{Float64}[]
    t_nnp1 = Complex{Float64}[]
    k_list = n_list.*(2*π/wavelength)
    kz_list = broadcast(*, k_list, broadcast(cos, θ_list))
    δ_list = kz_list.*d_list

    if polarization == "s"
        for i in 1:l-1
            push!(r_nnp1, r12_s(n_list[i], n_list[i+1], θ_list[i], θ_list[i+1]))
            push!(t_nnp1, t12_s(n_list[i], n_list[i+1], θ_list[i], θ_list[i+1]))
        end
    else
        for i in 1:l-1
            push!(r_nnp1, r12_p(n_list[i], n_list[i+1], θ_list[i], θ_list[i+1]))
            push!(t_nnp1, t12_p(n_list[i], n_list[i+1], θ_list[i], θ_list[i+1]))
        end
    end

    M_array = Matrix{Complex{Float64}}[]
    totalM = Matrix{Complex{Float64}}(I, 2, 2)
    for i in 1:l-1
        dn = δ_list[i]
        rn = r_nnp1[i]
        tn = t_nnp1[i]
        phase_mat = [exp(-1im * dn) 0; 0 exp(1im * dn)]
        rt_mat = [1/tn rn/tn; rn/tn 1/tn]
        push!(M_array, phase_mat * rt_mat)
        totalM = totalM * M_array[i]
    end

    r01 = r12_s(n_med, n_list[1], θi, θ_list[1])
    t01 = t12_s(n_med, n_list[1], θi, θ_list[1])
    M0 = [1/t01 r01/t01; r01/t01 1/t01]
    totalM = M0 * totalM
    return totalM
end

function plot_spectrum(wavelengths, polarization, θi, n_med, n_list, d_list, metal::Vector{Any})
    len = length(wavelengths)
    r = Complex{Float64}[]
    t = Complex{Float64}[]
    phase_r = Complex{Float64}[]

    for i in 1:len
        totalM = transfer_matrix(wavelengths[i], polarization, θi, n_med, n_list, d_list, metal)
        push!(t, 1/totalM[2,2])
        push!(r, totalM[2,1]/totalM[1,1])
    end

    R = broadcast(abs2, r)
    T = broadcast(abs2, t)

    return R, T
end

function plot_angle(wavelength, polarization, angles, n_med, n_list, d_list, metal::Vector{Any})
    len = length(angles)
    r = Complex{Float64}[]
    t = Complex{Float64}[]

    for i in 1:len
        totalM = transfer_matrix(wavelength, polarization, angles[i], n_med, n_list, d_list, metal)
        push!(t, 1/totalM[2,2])
        push!(r, totalM[2,1]/totalM[1,1])
    end

    R = broadcast(abs2, r)
    T = broadcast(abs2, t)

    return R, T
end

#%%
# Define parameters of multilayered material
θi = 0
θi = θi * π/180    # conversion to radian
n2, n1 = 1.45, 2.87
d2, d1 = 60, 150
N = 8
n_med = 1.0
c = 3.0 * 10^8
with_metal = [true, 168.26, 8935.2, 30, "left"]
without_metal = [false, 168.26, 8935.2, 30, "left"]

wavelengths = range(200, 800, step=1)

#%% Plotting
n_list, d_list = refractive_and_depths(N, n1, n2, d1, d2)
R, T = plot_spectrum(wavelengths, "s", θi, n_med, n_list, d_list, with_metal)
R0, T0 = plot_spectrum(wavelengths, "s", θi, n_med, n_list, d_list, without_metal)
R1, T1 = plot_spectrum(wavelengths, "s", θi, n_med, n_list, d_list, with_metal)
frequencies = broadcast(/, c, wavelengths)

#%%
plot(wavelengths, R0, label="R")
plot(wavelengths, T0, label="T")
xlabel("Wavelength (nm)")
#xlim([400, 700])
#plot(wavelengths, A0, label="A")
legend(loc="upper right")
title(L"$\theta_0 = 75\degree$")
#savefig("75.png")
#savefig("complex-phc-si.png")
#savefig("without_metal1.png")
#savefig("a-without-metal.png")

#%%
plot(wavelengths, R1)

#%% Absorption plot
function absorption(R, T)
    RT = broadcast(+, R, T)
    A = broadcast(-, 1, RT)
    return A
end

#%%
A0 = absorption(R0, T0)
A1 = absorption(R1, T1)
plot(wavelengths, A1, label="Au")
plot(wavelengths, A0, label="PhC")
legend()
#title("Absorption plot")
xlabel("Wavelength (nm)")
#savefig("au-absorption.png")

#%%
r0, t0 = plot_spectrum(wavelengths, "s", 0 * π/180, n_med, n_list, d_list, with_metal)
r1, t1 = plot_spectrum(wavelengths, "s", 30 * π/180, n_med, n_list, d_list, with_metal)
r2, t2 = plot_spectrum(wavelengths, "s", 45 * π/180, n_med, n_list, d_list, with_metal)
r3, t3 = plot_spectrum(wavelengths, "s", 60 * π/180, n_med, n_list, d_list, with_metal)
r4, t4 = plot_spectrum(wavelengths, "s", 75 * π/180, n_med, n_list, d_list, with_metal)
r5, t5 = plot_spectrum(wavelengths, "s", 85 * π/180, n_med, n_list, d_list, with_metal)

#%%
r0, t0 = plot_spectrum(wavelengths, "s", θi, n_med, n_list, d_list, [true, 168.26, 8935.2, 10, "left"])
r1, t1 = plot_spectrum(wavelengths, "s", θi, n_med, n_list, d_list, [true, 168.26, 8935.2, 20, "left"])
r2, t2 = plot_spectrum(wavelengths, "s", θi, n_med, n_list, d_list, [true, 168.26, 8935.2, 30, "left"])
r3, t3 = plot_spectrum(wavelengths, "s", θi, n_med, n_list, d_list, [true, 168.26, 8935.2, 40, "left"])
r4, t4 = plot_spectrum(wavelengths, "s", θi, n_med, n_list, d_list, [true, 168.26, 8935.2, 50, "left"])
r5, t5 = plot_spectrum(wavelengths, "s", θi, n_med, n_list, d_list, [true, 168.26, 8935.2, 60, "left"])

#%%
plot(wavelengths, r0, label="10 nm")
plot(wavelengths, r1, label="20 nm")
plot(wavelengths, r2, label="30 nm")
plot(wavelengths, r3, label="40 nm")
plot(wavelengths, r4, label="50 nm")
plot(wavelengths, r5, label="60 nm")
xlim([458, 590])
#ylim([-0.25, 1.2])
legend(loc="lower right")
xlabel("Wavelength (nm)")
ylabel("Reflectance")
#savefig("thickness-variation1.png")

#%%
plot(wavelengths, r0, label="0°")
plot(wavelengths, r1, label="30°")
plot(wavelengths, r2, label="45°")
plot(wavelengths, r3, label="60°")
plot(wavelengths, r4, label="75°")
plot(wavelengths, r5, label="85°")
xlim([480, 540])
#ylim([-0.25, 1.2])
legend()
xlabel("Wavelength (nm)")
ylabel("Reflectance")
#savefig("angle-variation.png")

#%%
plot(wavelengths, R1, label="R")
#xlim([500, 600])
plot(wavelengths, T1, label="T")
#plot(wavelengths, A1, label="A")
title(L"$\theta_0 = 0\degree$")
legend()
xlabel("Wavelength (nm)")
#savefig("a-with-metal.png")
#savefig("with-metal.png")


#%%
plot(wavelengths, R0, label="PhC")
plot(wavelengths, R1, label="Au")
legend()
xlabel("Wavelength (nm)")
ylabel("Reflectivity")
#savefig("high-n.png")
#xlim([400, 700])
#title("Tamm resonance")
#savefig("tamm.png")

#%%
plot(frequencies./1000, R, label="R")
plot(frequencies./1000, T, label="T")
legend()
xlabel("Frequency (THz)")

#%% Plotting at a fixed width of Au
θi = 40
θi = θi * π/180    # conversion to radian
n1, n2 = 1.45, 2.87
d1, d2 = 60, 150
N = 10
n_med = 1.0
c = 3.0 * 10^8
without_metal = [false, 168.26, 8935.2, 30, "left"]
with_metal = [true, 168.26, 8935.2, 30, "left"]
angles = range(0, 90, step=1)
θ = angles.*(π/180)    # conversion to radian

#%%
λ = 529
n_list, d_list = refractive_and_depths(N, n1, n2, d1, d2)
R1, T1 = plot_angle(λ, "s", angles, n_med, n_list, d_list, without_metal)
R2, T2 = plot_angle(λ, "s", angles, n_med, n_list, d_list, with_metal)

#%%
plot(angles, R1, label="PhC")
plot(angles, R2, label="Au coat")
legend()
xlabel("Incident angle")
ylabel("Reflectivity")


#%%
# For the case of Au on SiO2
θi = 40
θi = θi * π/180    # conversion to radian
n1, n2 = 1.45, 1.45
d1, d2 = 125*10^3, 125*10^3
N = 2
n_med = 1.0
c = 3.0 * 10^8
with_metal = [true, 168.26, 8935.2, 50, "left"]
without_metal = [false, 168.26, 8935.2, 30, "left"]

wavelengths = range(200, 800, step=1)

#%%
n_list, d_list = refractive_and_depths(N, n1, n2, d1, d2)
R0, T0 = plot_spectrum(wavelengths, "s", θi, n_med, n_list, d_list, without_metal)
R1, T1 = plot_spectrum(wavelengths, "s", θi, n_med, n_list, d_list, with_metal)

#%%
plot(wavelengths, R0, label="PhC")
plot(wavelengths, R1, label="Au")
legend()
xlabel("Wavelength (nm)")
ylabel("Reflectivity")

#%%
plot(wavelengths, R1, label="R")
plot(wavelengths, T1, label="T")
legend()
xlabel("Wavelength (nm)")


#%%
λ1 = 250
λ2 = 650
λ3 = 1800
angles = range(0, 720, step=1)
θ = angles.*(π/180)    # conversion to radian
R1, T1 = plot_angle(λ1, "s", θ, n_med, n_list, d_list, with_metal)
R2, T2 = plot_angle(λ2, "s", θ, n_med, n_list, d_list, with_metal)
R3, T3 = plot_angle(λ3, "s", θ, n_med, n_list, d_list, with_metal)

#%%
plot(angles, R1, label="250 nm")
plot(angles, R2, label="650 nm")
plot(angles, R3, label="1800 nm")
legend()
xlabel("Angle")
ylabel("Reflectivity")
