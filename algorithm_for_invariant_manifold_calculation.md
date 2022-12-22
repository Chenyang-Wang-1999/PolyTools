# Basic settings
- ODE equation:
$$\dot{x}=\mathcal{F}(x),x\in \mathscr{A}, \mathcal{F}\in\Gamma(T\mathscr{A})$$
- Invariant point $x_0$, near which $\mathcal{F}$ can be expanded to Taylor series of $u=x-x_0$. Denote the polynomial by $F$, then $F(0)=0$.
- Jacobi matrix $J$. Suppose $J$ can be block diagonalized as $$J=P\Lambda P^{-1}$$ where $$\Lambda=\begin{pmatrix} \Lambda_L & T \\0 & \Lambda_N \end{pmatrix}$$
- Solve the invariant manifold near $x_0$ is equivalent to find the functions $W:\Theta\rightarrow \mathscr{A}$ and $f\in\Gamma(T\Theta)$, so that $$F(W(s))=DW(s)f(s)$$

# Parameterization method
[1] À. Haro, M. Canadell, J.-L. Figueras, A. Luque, and J. M. Mondelo, The Parameterization Method for Invariant Manifolds: From Rigorous Results to Effective Computations, Vol. 195 (Springer International Publishing, Cham, 2016).

- Let $W(0)=0, \det(DW(0))\neq0$, then $f(0)=0$
- First order $W_1(s)=Ls$, $f_1(s)=\Lambda_Ls$, where $P=\begin{pmatrix}L & R\end{pmatrix}$.
- Iteration formula $$JW_k(s)-DW_k(s)\Lambda_Ls-Lf_k(s)=-E_k(s)$$where $$E_k(s)=\left[F(W_{<k}(s))\right]_k-\left[DW_{<k}(s)f_{<k}(s)\right]_k$$
- Transform by $P^{-1}$ $$P^{-1}JPP^{-1}W_k(s)-P^{-1}DW_k(s)\Lambda_Ls-P^{-1}Lf_k(s)=-P^{-1}E_k(s)$$ $\Rightarrow$ $$\begin{pmatrix}\Lambda_L & T \\ O & \Lambda_N\end{pmatrix}\xi_k - D\xi_k\Lambda_Ls - \begin{pmatrix}I_d\\O\end{pmatrix}f_k(s)=\eta_k(s)$$
- Tangent part and normal part: $$\xi_k=\begin{pmatrix} \xi^L_k \\ \xi_k^N \end{pmatrix}$$ then $$\begin{align*} \Lambda_L\xi^L_k+T\xi_k^N-D\xi_k^L\Lambda_Ls-f_k(s) &=\eta_k^L(s) \\ \Lambda_N\xi_k^N-D\xi_k^N\Lambda_Ls &=\eta^N_k(s) \end{align*}$$
- Terms of $\xi_k$. Suppose the monomial terms of $\xi_k$ is denoted as $\xi_m$, then $$D\xi_m\Lambda_Ls=\sum\limits_{j=1}^d m_j\lambda_j \xi_m :=m\lambda_L\xi_m$$  $$\begin{align*} \Lambda_L\xi^L_m+T\xi_m^N-m\lambda_L\xi_m^L-f_m(s) &=\eta_m^L(s) \\ \Lambda_N\xi_m^N-m\lambda_L\xi_m^N &=\eta^N_m(s) \end{align*}$$
- Elements of $\xi$ and $f$ $$\begin{align*} \lambda_j\xi^j_m+(T\xi_m^N)^j-m\lambda_L\xi_m^j-f_m^j(s) &=\eta_m^j(s),j\leq d \\ \Lambda_j\xi_m^j-m\lambda_L\xi_m^j &=\eta^j_m(s), j>d \end{align*}$$

# Code structure
- Input arguments: matrix $P$, matrix $T$, diagonal terms $\lambda_1,\dots, \lambda_n$, polynomial $F$
- Output arguments: polynomial $W$ and polynomial $f$
- Python interface
    - $P,T,\lambda$: use pybind11 Eigen interface with np.ndarray
    - $P^{-1}$: calculate internal by Eigen
    - Polynomials: add term by coefficients and orders, and utilities for python to read the terms. In python, convert the coefficient and order data to sympy expression
    - Evaluation for numpy array
- Algorithm:
    - Initialize all data, add $f_1$ and $W_1$
    - For $k=2$ to $k=k_{\max}-1$
        - Calculate $E_k$ (Homogen) 
        - Calculate $\eta_k$
        - Update $\xi_k^N$ (for each term in $\eta_k$, add the terms to $\xi_k$.)
        - Update $T\xi_k^N$ if $T$ is non empty
        - Update $\xi^L_k$ and $f_k$
        - Multiply $P$ on the left of $\xi_k^L$, destructively add them to $W$