# Milky Way Potentials

Following @borukhovetskaya+2022, we define a multi-component MW potential as follows (in code units: 1 kpc, 10^10 Msun)



| Component  | Values                    |      |
| ---------- | ------------------------- | ---- |
| thin disk  | M=5.9, a=3.9, b=0.31      |      |
| thick disk | M=2, a=4.4, b=0.92        |      |
| bulge      | M = 2.1, a = 1.3          |      |
| halo       | Mvir=115, r=20.2, c=9.545 |      |





instead the @vasily+2021 potential is

| Component         | Values                                      |      |
| ----------------- | ------------------------------------------- | ---- |
| spherical bulge   | $\Gamma=1.8$, $r_b=0.2$, $u_b=1.8$, $M=1.2$ |      |
| disk (isothermal) | M=5, r_d=3, h_d=0.4                         |      |
| halo              |                                             |      |

## 



## Notes about derivations

For my sake, here are the equations to derive each of the following relations. 

Given a density profile, the contained mass is $M(r) = \int_0^r \rho dV$. The acceleration is $a = -\nabla \Phi$. 



## NFW (Halo)

### Aside:Defined quantities

The NFW is parameterized in terms of 

with $c=r_{200}/r_s$

The NFW halo is written in terms of $M_{200}$, $R_{200}$. These are the radius and mass of a sphere with 200 times the average critical density, $\rho_{crit} = 3H^2/8\pi G=127.35M_\odot/{\rm kpc}^3 $, so
$$
\rho_{200} = 200\rho_{crit} = \frac{M_{200}}{4/3\ \pi R_{200}^3}
$$

$$
R_{200} = \sqrt[3]{\frac{1}{200}\frac{3M_{200}}{4\pi \rho_{\rm crit}}}
$$


$$
M_{200} = \frac{4\pi}{3} R_{200}^3\ \rho_{200}
$$


nother quantity which turns up often here is
$$
A(x) \equiv \log (1+x) - \frac{x}{1+x}
$$


We will also just define now
$$
x \equiv r/r_s
$$
as all quantities are simply scaled by the scale radius.

For example, Asya uses $M_{200} = 1.04e10$ and $c=12.5$ for Fornax, giving $M_s = 0.119\times10^{10}\,{\rm M}_\odot$ and $R_{s}=3.68\,$kpc.

### Density

From the @nfw1996 paper, 
$$
\rho(x) =  \frac{\rho_s/3}{x\ (1+x)^2}
$$
where
$$
\rho_s = \frac{c^3}{A(c)} \rho_{200}
$$
and $A(c)$ is as above. The characteristic density can also be written in terms of scale mass,  $M_s = M_{200}/{A(c)}$  (see below), giving
$$
\rho_s =  \frac{M_s}{4\ \pi\, {r_s}^3}
$$
### Contained mass

By integrating $\rho$, the contained mass interior to $r$ is
$$
M(x) = M_s\ A(x)
$$
where
$$
M_s \equiv 4\pi/3 \rho_s r_s^3 = \frac{M_{200}}{A(c)}
$$
 so $M(x) = M_s  A(x)$. Note that this differs from the mass contained within $r_s$ by $A(1) \approx 0.193147$.  

### Potential

From Poisson's equation, 
$$
\Phi(x) = -\frac{G M_s}{r_s} \ \frac{1}{x} \ \ln\left(1+x\right)
$$

 (The constant of integration are chosen such that $\Phi(0) = \Phi_0$ and $\Phi$ vanishes at $\infty$.)

### Acceleration 

By differentiating the potential,
$$
\vec{a}(x) = -\frac{G M_s}{r_s^2} \frac{A(x)}{x^2} \hat{r}
$$
Also is consistent with using contained mass ($F=G M(r)/r^2$)

### Circular velocity

if $x = r/r_s$, then 
$$
\left(V_c/V_{200}\right)^2 = \frac{A(x)/x}{A(c)/c}
$$

from the radius and maximum circular velocity, 
$$
R_{\rm circ}^{\rm max} = \alpha\ r_s
$$

where $\alpha\approx2.16258$



## Truncated NFW



## Hernquist (Bulge)

The @hernquist1990 density profile for the galactic bulge is parameterized in terms of a characteristic mass, $M$, and radius $a$. 

Density Profile
$$
\rho(r) = \frac{M}{2\pi} \frac{a}{r} \frac{1}{(r+a)^3}
$$
Mass profile
$$
M(r) = M \frac{r^2}{(r+a)^2}
$$
Potential:
$$
\Phi(r) = - \frac{GM}{r + a}
$$
The acceleration is then
$$
\vec{a} = - \hat{r} \frac{G\,M}{\left(r + a\right)^2}
$$

Circular velocity

## Miyamoto-Nagai (Bulge)

$$
\Phi(r) = - \frac{M}{\sqrt{R^2 + b^2}}
$$

(special case of $a=0$ for Miyamoto Nagai Disk)

## Miyamoto-Nagai (Disk)

for the MW disk (thin and thick), the potential is derived from @miyamoto1975, 
$$
\Phi(R, z) = \frac{-GM}{\left(R^2 + \left[a + \sqrt{z^2 + b^2}\right]^{2}\right)^{1/2}}
$$
for convenience, we define the helper variables, $S = \sqrt{z^2 + b^2}$ and $D = \sqrt{R^2 + (a+S)^2}$, so the potential simplifies to 
$$
\Phi(R, z) = \frac{-GM}{D}
$$

$$
a_{\rm disk} = -\nabla \Phi = -\left(\hat{r}\frac{\partial}{\partial r} + \hat{z} \frac{\partial}{\partial{z}}\right) \Phi
$$

$$
a_{\rm disk} = -\frac{G\,M}{D^3}\,R\,\hat{R} - \frac{G\,M}{D^3} \left(\frac{a}{S}+1\right) \,z\,\hat{z}
$$



density
$$
\rho = \frac{b^2 M}{4\pi} \frac{aR^2 + (a+3r)(a + r)^2}{S^3 D^5}
$$




## Powerlaw Bulge (Vasily)


$$
\rho	\propto (1 + r/r_b)^{-\Gamma} \exp(-(r/u_b)^2)
$$



## Logarithmic

$$
\Phi = v^2/2 \ln\left(x^2 + y^2 + (z/q)^2 + d^2\right)
$$

acceleration
$$
\vec{a} = -\frac{v^2}{x^2 + y^2 + (z/q)^2 + d^2} 
\left(
	\vec{x} + \vec{y} + \vec{z} / q^2
\right)
$$




## Isothermal (disk)

$$
\rho	\sim \exp(-R/R_d)\ {\rm sech}^2(z/(2h))
$$



## $\alpha\beta\gamma$ profile (bulge)

 @zhau1996.


$$
\rho = r^{-\gamma} (1 + r^{\alpha})^{(\beta - \gamma)/\alpha}
$$


potential is compl

## $\alpha\beta\gamma$​​ profile (bulge?)


$$
\rho_h = (s / r_h)^{-\gamma} (1 + (s/r_h)^\alpha)^{(\gamma - \beta)/\alpha} \exp(-(s/u_h)^\eta)\\
s \equiv (pq)^{1/3} \sqrt{x^2 + (y/p)^2 + (z/q)^2}
$$

where $s$ is sphericalized (or hyperelliptical) radius 

with exponential drop off $u_h = 200$, $\eta=2$. Becomes as published with $u_h \to \infty$





# 2D profiles

## Exponential

$$
\Sigma(R) = \Sigma_0 \exp(-R/R_s)
$$


$$
M(R) =  \pi  R_s^2 \Sigma_0 \left(1 - x \exp(-x) - \exp(-x)\right)
$$


Solving $M(R) = 1/2M_{\rm tot}$,  the half mass radius is $1.6783\,R_s$​. 

for the 3D exponential, 
$$
M(r) = 4\pi r_s^3 \rho_0 \left(2 - 2e^{-x} - 2xe^{-x} - x^2\,e^{-x}\right)
$$
which gives $r_h = 2.674 r_s$. So, given $R_s$, $r_s = 0.62762\,R_s$. (Verify this!)
