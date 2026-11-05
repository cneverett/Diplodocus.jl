# Synchrotron Radiation
The synchrotron module allows for the evaluation of the rate of change of photon population (``\partial_t f(p_1,\cos\theta_1)``) due to synchrotron emissions from a distribution (``f(p_2,\cos\theta_2)``) charged particles, with mass ``m_2`` and charge ``q_2``, in a uniform magnetic field of strength ``B`` aligned with the axis of symmetry.

The rate of change of photon distribution function is given by:
```math
\left(\frac{\partial f(\boldsymbol{p}_1)}{\partial t}\right)_{Sync}=\frac{Z}{B}\frac{3c^5\sigma_Tm_e^5}{4\pi\hbar^3\mu_0e^3}\int\mathrm{d}^3\boldsymbol{p}_2\frac{1}{p_1p^0_2}\left[\left(\frac{p^0_2\cos\theta_1-p_2\cos\theta_2}{\sin\theta_1}\right)^2J_{\omega/\omega_0}(x)^2+(p_2\sin\theta_2)^2 J'_{\omega/\omega_0}(x)^2\right]f(\boldsymbol{p}_2).
```
where ``J`` are the Bessel functions (``J'`` refer to the derivative with respect to the argument), ``\frac{\omega}{\omega_0}=\frac{m_e^2c^2p_1}{Z\hbar eB}\left(p^0_2-p_1\cos\theta_1\cos\theta_2\right)`` and ``x=\frac{\omega}{\omega_0}=\frac{m_e^2c^2p_1}{Z\hbar eB}p_1p_2\sin\theta_1\sin\theta_2``. ``Z`` and ``\mu`` are the normalised charge and mass of the emitting particle, normalised by the fundamental charge and electron mass respectively, all other constants have their standard (SI) values.

Making the standard set of assumptions laid out in [BoltzmannCollisionIntegral.jl](@ref) this can be re-written in the form
```math
\left(\frac{\partial f_{1,ij}}{\partial t}\right)_{sync}=S_{ijkl}f_{2,kl}
```
with 
```math
S_{ijkl}=\frac{2\pi}{\Delta p_i\Delta u_j}\int_{\Delta p_i\Delta u_j}p_1^2\mathrm{d}p_1\mathrm{d}u_1 \int \mathrm{d}p_2\mathrm{d}u_2~\left|\frac{Z}{B}\frac{3c^5\sigma_Tm_e^5}{4\pi\hbar^3\mu_0e^3}\right|\frac{1}{p_1p^0_2}\left[\left(\frac{p^0_2 u_1-p_2u_2}{\sqrt{1-u^2_1}}\right)^2J_{\omega/\omega_0}(x)^2+p_1^2(1-u_2^2) J'_{\omega/\omega_0}(x)^2\right]
```
where ``u=\cos\theta``.

!!! note
    Evaluation is normalised by ``c\sigma_T`` to match the normalisation of the binary collision integral.

The values of ``S_{ijkl}`` are then evaluated using a Monte Carlo method. To improve accuracy when ``\omega/\omega_0\ge 10^3``, the following approximation of the Bessel functions are made: 
```math
J_n(x)\approx\frac{1}{\pi\sqrt{3}}\left(1-\frac{x^2}{n^2}\right)^{1/2}K_{1/3}\left(\frac{n}{3}\left(1-\frac{x^2}{n^2}\right)^{3/2}\right),
```
```math
J'_n(x)\approx\frac{1}{\pi\sqrt{3}}\left(1-\frac{x^2}{n^2}\right)K_{2/3}\left(\frac{n}{3}\left(1-\frac{x^2}{n^2}\right)^{3/2}\right).
```