module NummatDn03

export theta_at_time, find_period, period_vs_energy, harmonic_theta, demo_main
using Printf

# helper funkcije
# max abs vrednost vektorja (inf norm)
function infnorm(v::AbstractVector{<:Real})
    m = 0.0
    @inbounds for i in eachindex(v)
        x = abs(float(v[i]))
        if x > m
            m = x
        end
    end
    return m
end

# napaka za adaptive step size
function scaled_error(y5::AbstractVector{<:Real}, y4::AbstractVector{<:Real}, atol::Real, rtol::Real)
    err = 0.0
    @inbounds for i in eachindex(y5)
        sc = atol + rtol*max(abs(float(y5[i])), abs(float(y4[i])))   # scale factor
        e  = abs(float(y5[i] - y4[i])) / sc
        if e > err
            err = e
        end
    end
    return err
end

#DOPRI5 integrator
#Butcherjev table

mutable struct Dopri5Result
    t::Vector{Float64}
    y::Vector{Vector{Float64}}
end

#en step DOPRI5, vrača oceno 5. reda, 4. reda in k7 (za FSAL)
function dopri5_step!(f, t, y, h, k1_cache=nothing)

    # koeficienti
    a21=1//5
    a31=3//40; a32=9//40
    a41=44//45; a42=-56//15; a43=32//9
    a51=19372//6561; a52=-25360//2187; a53=64448//6561; a54=-212//729
    a61=9017//3168; a62=-355//33; a63=46732//5247; a64=49//176; a65=-5103//18656

    a71=35//384; a72=0//1; a73=500//1113; a74=125//192; a75=-2187//6784; a76=11//84
    b1=a71; b2=a72; b3=a73; b4=a74; b5=a75; b6=a76; b7=0//1

    # 4th order weights
    bs1=5179//57600; bs2=0//1; bs3=7571//16695; bs4=393//640; bs5=-92097//339200; bs6=187//2100; bs7=1//40

    c2=1//5; c3=3//10; c4=4//5; c5=8//9; c6=1//1; c7=1//1

    # k1 je lahko podan iz prejšnjega stepa
    k1 = k1_cache === nothing ? f(t,y) : k1_cache

    # naprej po korakih
    y2 = @. y + h*Float64(a21)*k1
    k2 = f(t + Float64(c2)*h, y2)

    y3 = @. y + h*(Float64(a31)*k1 + Float64(a32)*k2)
    k3 = f(t + Float64(c3)*h, y3)

    y4 = @. y + h*(Float64(a41)*k1 + Float64(a42)*k2 + Float64(a43)*k3)
    k4 = f(t + Float64(c4)*h, y4)

    y5 = @. y + h*(Float64(a51)*k1 + Float64(a52)*k2 + Float64(a53)*k3 + Float64(a54)*k4)
    k5 = f(t + Float64(c5)*h, y5)

    y6 = @. y + h*(Float64(a61)*k1 + Float64(a62)*k2 + Float64(a63)*k3 + Float64(a64)*k4 + Float64(a65)*k5)
    k6 = f(t + Float64(c6)*h, y6)

    # ocena 5. reda
    y5th = @. y + h*(Float64(b1)*k1 + Float64(b2)*k2 + Float64(b3)*k3 + Float64(b4)*k4 + Float64(b5)*k5 + Float64(b6)*k6)
    k7   = f(t + Float64(c7)*h, y5th)

    # ocena 4. reda
    y4th = @. y + h*(Float64(bs1)*k1 + Float64(bs2)*k2 + Float64(bs3)*k3 + Float64(bs4)*k4 + Float64(bs5)*k5 + Float64(bs6)*k6 + Float64(bs7)*k7)

    return y5th, y4th, k7
end

#glavna funkcija dopri5 (adaptive)
function dopri5!(f, t0, y0, tf; rtol=1e-9, atol=1e-12, h_init=nothing, h_min=1e-12, h_max=Inf, maxsteps=10^7)
    t=float(t0); tf=float(tf); y=Float64.(y0)
    h = h_init === nothing ? min(0.1, abs(tf-t)/10) : float(h_init)

    Ts = [t]; Ys = [copy(y)]
    k1c = nothing; steps=0

    while (t < tf && h > 0) || (t > tf && h < 0)
        steps+=1
        if steps>maxsteps; error("preveč korakov..."); end
        if (t+h>tf && h>0)||(t+h<tf && h<0); h=tf-t; end

        y5, y4, k7 = dopri5_step!(f,t,y,h,k1c)
        err = scaled_error(y5,y4,atol,rtol)

        if err <= 1.0 || abs(h)<=h_min
            # sprejmemo step
            t += h; y=y5
            push!(Ts,t); push!(Ys,y)
            k1c=k7
            fac = err==0.0 ? 5.0 : clamp(0.9*err^(-0.2),0.2,5.0)
            h = sign(h)*min(abs(h)*fac,h_max)
        else
            fac = clamp(0.9*err^(-0.2),0.1,0.5)
            h = sign(h)*max(abs(h)*fac,h_min)
            k1c=nothing
        end

        if t==tf; break; end
    end
    Dopri5Result(Ts,Ys)
end

# solve samo do enega časa
function solve_to(f, t0, y0, t1; rtol=1e-9, atol=1e-12, h_init=nothing)
    if t0==t1
        return Dopri5Result([float(t0)],[Float64.(y0)])
    end
    dopri5!(f,t0,y0,t1;rtol=rtol,atol=atol,h_init=h_init)
end

#nihalo
struct PendulumParams; g::Float64; l::Float64; end

# right hand side za pendulum
function pendulum_rhs(t,y,p::PendulumParams)
    θ=y[1]; ω=y[2]
    [ω, -(p.g/p.l)*sin(θ)]
end



"""
theta_at_time(t; θ0, ω0, g=9.81, l=1.0)

Vrne odmik nihala θ(t) pri času `t`.
Parametri: začetni odmik `θ0` (rad), začetna hitrost `ω0` (rad/s), `g` in `l`.
Uporablja lastno DOPRI5 integracijo.
"""
# θ(t)
function theta_at_time(t; θ0,ω0,g=9.81,l=1.0,rtol=1e-9,atol=1e-12)
    p=PendulumParams(float(g),float(l))
    f=(τ,y)->pendulum_rhs(τ,y,p)
    sol=solve_to(f,0.0,[float(θ0),float(ω0)],float(t);rtol=rtol,atol=atol)
    sol.y[end][1]
end

energy(θ,ω;g=9.81,l=1.0)= g*l*(1-cos(θ))+0.5*(l^2)*ω^2




"""
harmonic_theta(t, θ0, ω0; g=9.81, l=1.0)

Linearni (majhni odmiki) model nihala: θ(t) ≈ θ0 cos(√(g/l) t) + (ω0/√(g/l)) sin(√(g/l) t).
Služi kot primerjava, ni “točen” za velike amplitude.
"""
function harmonic_theta(t,θ0,ω0;g=9.81,l=1.0)
    wn=sqrt(g/l)
    θ0*cos(wn*t)+(ω0/wn)*sin(wn*t)
end

advance_to(f,tL,yL,tR;rtol=1e-9,atol=1e-12)=
    solve_to(f,tL,yL,tR;rtol=rtol,atol=atol).y[end]

# refine zero crossing za θ=0
function refine_zero_theta(f,tL,yL,tR,yR;rtol=1e-9,atol=1e-12,tol_time=1e-12)
    θL=yL[1]; θR=yR[1]
    while abs(tR-tL)>tol_time
        tM=0.5*(tL+tR)
        yM=advance_to(f,tL,yL,tM;rtol=rtol,atol=atol)
        θM=yM[1]
        if θL==0.0; return tL,yL
        elseif θR==0.0; return tR,yR
        end
        if sign(θL)!=sign(θM); tR,yR,θR=tM,yM,θM
        else; tL,yL,θL=tM,yM,θM; end
    end
    tZ=0.5*(tL+tR)
    yZ=advance_to(f,tL,yL,tZ;rtol=rtol,atol=atol)
    return tZ,yZ
end



"""
find_period(θ0, ω0; g=9.81, l=1.0)

Oceni periodo nihanja z iskanjem dveh zaporednih prehodov skozi θ=0 (smer navzgor, ω>0).
Vrne čas v sekundah. Za večje amplitude je T > 2π√(l/g).
"""
# najdi periodo prek preseka θ=0, samo naraščajoč prehod
function find_period(θ0,ω0;g=9.81,l=1.0,rtol=1e-9,atol=1e-12,maxT=1e4)
    p=PendulumParams(float(g),float(l))
    f=(τ,y)->pendulum_rhs(τ,y,p)

    t=0.0; y=[float(θ0),float(ω0)]
    h=0.05; crossings=Float64[]; last_t=t; last_y=copy(y)

    while t<maxT && length(crossings)<2
        sol=solve_to(f,t,y,t+h;rtol=rtol,atol=atol)
        t=sol.t[end]; y=sol.y[end]

        θL=last_y[1]; θR=y[1]
        if (θL<0.0)&&(θR>0.0)
            tZ,yZ=refine_zero_theta(f,last_t,last_y,t,y;rtol=rtol,atol=atol,tol_time=1e-12)
            if yZ[2]>0.0; push!(crossings,tZ); end
            t,y=tZ,yZ
            sol2=solve_to(f,t,y,t+1e-6;rtol=rtol,atol=atol) # micro step
            t=sol2.t[end]; y=sol2.y[end]
        end
        last_t=t; last_y=copy(y)
        h=min(0.2,max(0.01,h))
    end

    if length(crossings)<2; error("no 2 crossings :("); end
    crossings[end]-crossings[end-1]
end


"""
period_vs_energy(θlist; g=9.81, l=1.0)

Za seznam začetnih odmikov `θlist` (ω0=0) vrne (E, T), kjer je
E = g*l*(1 - cos θ0), T = perioda za ta θ0.
Uporabno za graf T(E). 
"""
# loop čez več θ0 in vrni T(E)
function period_vs_energy(θlist;g=9.81,l=1.0,rtol=1e-9,atol=1e-12)
    Es=Float64[]; Ts=Float64[]
    for θ0 in θlist
        E=energy(float(θ0),0.0;g=g,l=l)
        T=find_period(θ0,0.0;g=g,l=l,rtol=rtol,atol=atol)
        push!(Es,E); push!(Ts,T)
    end
    return Es,Ts
end

#demo
function demo_main()
    g=9.81; l=1.0

    println("primer 1: theta(t) pri θ0=0.5, t=1s")
    θt=theta_at_time(1.0;θ0=0.5,ω0=0.0,g=g,l=l)
    println("θ(1) ≈ ",@sprintf("%.12f",θt))

    println("\nprimer 2: perioda za θ0=1 rad")
    T=find_period(1.0,0.0;g=g,l=l)
    println("T ≈ ",@sprintf("%.12f s",T))

    println("\nprimer 3: primerjava s harmonicnim")
    θ_num=theta_at_time(2.0;θ0=0.1,ω0=0.0,g=g,l=l)
    θ_har=harmonic_theta(2.0,0.1,0.0;g=g,l=l)
    println("θ_num(2s) ≈ ",@sprintf("%.12f",θ_num),", θ_har(2s) ≈ ",@sprintf("%.12f",θ_har))

    println("\nprimer 4: T(E) za θ0 v [0.05,2.9]")
    θs=collect(range(0.05,stop=2.9,length=20))
    Es,Ts=period_vs_energy(θs;g=g,l=l)
    for (E,T) in zip(Es,Ts)
        println(@sprintf("E=%.6f, T=%.6f",E,T))
    end
    println("\nkonc")
end

end # module
