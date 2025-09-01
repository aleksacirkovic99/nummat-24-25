using Test
include(joinpath(@__DIR__, "..", "src", "NummatDn03.jl"))
using .NummatDn03

@testset "Matematično nihalo (DOPRI5)" begin
    g = 9.81; l = 1.0

    # 1) majhen odmik ≈ harmonično nihalo
    θ0 = 0.05; ω0 = 0.0; t = 3.0
    θ_num = theta_at_time(t; θ0=θ0, ω0=ω0, g=g, l=l)
    θ_lin = harmonic_theta(t, θ0, ω0; g=g, l=l)
    @test isapprox(θ_num, θ_lin; rtol=1e-4)

    # 2) perioda blizu teoretične (s popravkom za majhno amplitudo)
    T_num = find_period(θ0, ω0; g=g, l=l)
    T0 = 2π*sqrt(l/g)            # linearna perioda
    T_corr = T0*(1 + θ0^2/16)    # prva korekcija za majhen odmik
    @test isapprox(T_num, T_corr; rtol=1e-5)

    # 3) večja amplituda → daljša perioda
    T1 = find_period(0.2, 0.0; g=g, l=l)
    T2 = find_period(1.0, 0.0; g=g, l=l)
    @test T2 > T1
end
