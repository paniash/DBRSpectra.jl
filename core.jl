using LinearAlgebra
using PyPlot
pyplot()

#%%
# Define parameters of multilayered material
θi = 0
n1, n2 = 1.5, 2.5
d1, d2 = 50, 90
N = 10
n_med = 1.0


#%%
wavelengths = range(550, 1210, step=1)

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
    θ_list = Float64[]
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

function M_matrix(wavelength, polarization, θi, n_med, n_list, d_list)
    push!(n_list, 1)
    push!(d_list, 10)

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

function plot_spectrum(wavelengths, polarization, n_med, n_list, d_list)
    len = length(wavelengths)
    r = Complex{Float64}[]
    t = Complex{Float64}[]
    phase_r = Complex{Float64}[]

    for i in 1:len
        totalM = M_matrix(wavelengths[i], polarization, θi, n_med, n_list, d_list)
        push!(t, 1/totalM[2,2])
        push!(r, totalM[2,1]/totalM[1,1])
        push!(phase_r, phase_of_reflectivity(totalM[2,1]/totalM[1,1]))
    end

    R = broadcast(abs2, r)
    T = broadcast(abs2, t)

    return wavelengths, phase_r, R, T
end

#%% Plotting
n_list, d_list = refractive_and_depths(N, n1, n2, d1, d2)
wavelengths, phase_r, R, T = plot_spectrum(wavelengths, "s", n_med, n_list, d_list)

#%%
plot(wavelengths, R, label="R")
plot(wavelengths, T, label="T")
legend()
xlabel("Wavelength (nm)")
