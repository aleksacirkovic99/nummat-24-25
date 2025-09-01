# DN03 – Matematično nihalo
**Avtor:** Aleksa Ćirković (63230486)
**Predmet:** Numerična matematika  
**Jezik:** Julia 1.11

## 1. Opis naloge
Računamo rešitev nedušenega matematičnega nihala:
\[
\ddot{\theta}(t) + \frac{g}{l}\sin\theta(t) = 0,\quad \theta(0)=\theta_0,\ \dot{\theta}(0)=\omega_0.
\]

Zahteva: implementacija adaptivne metode DOPRI5, funkcija za \(\theta(t)\), numerična perioda in odvisnost T(E).

## 2. Opis rešitve
Sistem pretvorimo v sistem 1. reda:  
\(\dot{\theta} = \omega,\ \dot{\omega} = -(g/l)\sin\theta.\)  
Za integracijo uporabimo metodo Dormand–Prince 5(4)7 (lastna implementacija, brez ODE paketov).  
Perioda se izračuna s pomočjo preseka \(\theta=0\), ko je \(\omega>0\).

## 3. Rezultati
Primeri, ki jih da `demo.jl`:
- \(\theta(1s)\) za \(\theta_0=0.5\).  
- Perioda pri \(\theta_0=1.0\).  
- Primerjava z linearnim modelom pri majhnih odmikh.  
- T(E) tabela za različne amplitude.

*(Sem prilepiš screenshote iz terminala ali graf iz demota — “slikice prosim”)*

## 4. Matematični komentar
Za majhne odmike:
\[
T \approx 2\pi\sqrt{l/g}\Big(1+\frac{\theta_0^2}{16}\Big).
\]
Numerični rezultati se skladajo s teorijo: večja amplituda → daljša perioda.

## 5. Navodila za uporabo
```julia
] activate .
include("demo.jl")

## 6. Testi
] activate .
] test
