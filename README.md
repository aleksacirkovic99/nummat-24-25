# Aleksa Ćirković (63230486)


# DN03 — Matematično nihalo (Julia, lastni DOPRI5)
Implementacija nedušenega matematičnega nihala z **lastno** adaptivno metodo **DOPRI5 (Dormand–Prince 5(4)7)** brez uporabe ODE paketov.  
Zagotovljene funkcije:
- `theta_at_time(t; θ0, ω0, g, l)` – odmik nihala ob času `t`.
- `find_period(θ0, ω0; g, l)` – numerična perioda prek Poincaréjevega preseka `θ=0` z `ω>0`.
- `period_vs_energy(θ0_list; g, l)` – podatki za T(E).
- `harmonic_theta(t, θ0, ω0; g, l)` – analitični linearni model za primerjavo.

# Zahteve
- Julia **1.11** (testirano na 1.11.6)
- Standardne knjižnice: `Printf`, `Test` (dodani v projekt)

# Struktura
NummatDn03/
─ Project.toml
─ Manifest.toml
─ README.md
─ src/
 ─ NummatDn03.jl # modul z DOPRI5 in API funkcijami
─ demo.jl # primer uporabe (kliče demo_main())
─ test/
 ─ runtests.jl # enotski testi


# Hiter zagon (demo)
V Julia REPL-ju v korenu projekta:
```julia
] activate .
include("demo.jl")   # pokliče demo_main() in izpiše primere

# Test
] activate .
] test



# Kratek opis
theta_at_time(t; θ0, ω0, g=9.81, l=1.0)  :: Float64
find_period(θ0, ω0; g=9.81, l=1.0)       :: Float64
period_vs_energy(θ0_list; g=9.81, l=1.0) :: Tuple{Vector{Float64},Vector{Float64}}
harmonic_theta(t, θ0, ω0; g=9.81, l=1.0) :: Float64
