Mathematical Foundations, Theoretical Mechanics, and Numerical Stability Boundaries of Computational Algorithms
===============================================================================================================

## Nonlinear Root-Finding and Functional Operator Inversion

### Mathematical Mechanics of Newton-Raphson and Newton-Kantorovich Frameworks

The classical Newton-Raphson method addresses the fundamental problem of identifying a zero $x^* \in \mathbb{R}^n$ of a continuously differentiable mapping $F: \mathbb{R}^n \to \mathbb{R}^n$ (Ortega & Rheinboldt, 1970; Polyak, 2007). The foundation of this technique rests on a local first-order Taylor expansion about a given iterate $x_k$:

$$F(x) = F(x_k) + F'(x_k)(x - x_k) + \mathcal{R}_1(x)$$

Neglecting the integral remainder term $\mathcal{R}_1(x) = \int_0^1 (1-t) F''(x_k + t(x - x_k))(x - x_k)^2 \, dt$ yields an affine model of the operator, denoted $M_k(x) = F(x_k) + F'(x_k)(x - x_k)$ (Argyros, 1998; Kantorovich, 1948). Enforcing the condition $M_k(x_{k+1}) = 0$ yields the iterative scheme:

$$x_{k+1} = x_k - [F'(x_k)]^{-1} F(x_k)$$

This iteration requires the resolution of the linear system $F'(x_k) s_k = -F(x_k)$ followed by the update $x_{k+1} = x_k + s_k$ at every discrete stage (Polyak, 2007).

In infinite-dimensional functional analysis, this concept generalizes to the Newton-Kantorovich framework, which solves operator equations $F(x) = 0$ defined on open subsets $\Omega$ of a Banach space $X$ with values in a Banach space $Y$ (Ciarlet & Mardare, 2014; Kantorovich & Akilov, 1982). Originating from Leonid Kantorovich's 1948 contributions, this formulation avoids the circular assumption that an exact solution $x^*$ exists prior to convergence analysis (Kantorovich, 1948; Polyak, 2007). Instead, semi-local convergence conditions rely exclusively on structural quantities measured at the initial point $x_0 \in \Omega$ (Argyros, 1998; Ferreira & Svaiter, 2012).

Theorem (Newton-Kantorovich):

Let $X$ and $Y$ be Banach spaces, $\Omega \subset X$ an open convex set, and $F: \Omega \to Y$ a Fréchet differentiable mapping (Ciarlet & Mardare, 2014). Assume that at an initial vector $x_0 \in \Omega$, the linear operator $\Gamma_0 = [F'(x_0)]^{-1} \in \mathcal{L}(Y, X)$ exists as a bounded inverse, satisfying the following bounds (Ciarlet & Mardare, 2014; Magreñán, 2016):

1. $\|\Gamma_0\|_{\mathcal{L}(Y, X)} \le \beta$

2. $\|\Gamma_0 F(x_0)\|_X \le \eta$

3. $\|F'(x) - F'(y)\|_{\mathcal{L}(X, Y)} \le L \|x - y\|_X, \quad \forall x, y \in \Omega$

Define the dimensionless Kantorovich parameter:

$$h = \beta \eta L$$

If $h \le \frac{1}{2}$ and the closed ball $\bar{B}(x_0, r_-) \subseteq \Omega$, where the radius is defined by the minor root of the scalar quadratic majorant $p(t) = \frac{1}{2} \beta L t^2 - t + \eta$ (Argyros, 1998; Kantorovich, 1948):

$$r_- = \frac{1 - \sqrt{1 - 2h}}{\beta L} = \frac{2\eta}{1 + \sqrt{1 - 2h}}$$

then the sequence of Newton iterates $\{x_k\}$ is well-defined, remains entirely within the closed domain $\bar{B}(x_0, r_-)$, and converges strongly to a unique solution $x^* \in \bar{B}(x_0, r_-)$ (Ciarlet & Mardare, 2014; Kantorovich & Akilov, 1982). Furthermore, this solution is geometrically unique within the larger open ball $B(x_0, r_+)$, where $r_+ = \frac{1 + \sqrt{1 - 2h}}{\beta L}$, collapsing to uniqueness in $\bar{B}(x_0, r_-)$ precisely at the critical threshold $h = \frac{1}{2}$ (Ciarlet & Mardare, 2014).

The mechanism underlying this proof depends on constructing the scalar majorizing sequence $t_{k+1} = t_k - \frac{p(t_k)}{p'(t_k)}$, initialized at $t_0 = 0$ (Argyros, 1998; Ferreira & Svaiter, 2012). By invoking the Banach Lemma on invertible operators, whenever $\|I - \Gamma_0 F'(x)\| < 1$, the operator $F'(x)$ remains invertible with its norm controlled by the scalar derivative $|p'(t)|^{-1}$ (Argyros, 1998; Kantorovich & Akilov, 1982). An induction argument confirms the componentwise domination:

$$\|x_{k+1} - x_k\| \le t_{k+1} - t_k$$

Because $p(t)$ is convex, the scalar iterates converge monotonically to $r_-$, forcing $\{x_k\}$ to form a Cauchy sequence in $X$ that converges to $x^*$ at an asymptotic quadratic rate (Argyros, 1998; Ferreira & Svaiter, 2012):

$$\|x_{k+1} - x^*\| \le \frac{\beta L}{2(1 - \beta L r_-)} \|x_k - x^*\|^2$$

### Brent's Hybrid Root-Finding Architecture

The van Wijngaarden-Dekker-Brent algorithm provides root containment for continuous scalar functions $f: [a, b] \to \mathbb{R}$ where $f(a)f(b) \le 0$ without requiring evaluation of the derivative $f'(x)$ (Brent, 1973; Dekker, 1969). The method couples the linear reliability of bisection with the superlinear acceleration of open interpolation schemes (Brent, 1973).

The primary acceleration step is Inverse Quadratic Interpolation (IQI) (Brent, 1973). Given three distinct historical support coordinates $(a_k, f(a_k))$, $(b_k, f(b_k))$, and $(c_k, f(c_k))$, where $b_k$ denotes the current best approximation such that $|f(b_k)| \le |f(a_k)|$, the inverse relationship $x = f^{-1}(y)$ is modeled via a second-order Lagrange polynomial:

$$x = \frac{(y - f(b_k))(y - f(c_k))}{(f(a_k) - f(b_k))(f(a_k) - f(c_k))} a_k + \frac{(y - f(a_k))(y - f(c_k))}{(f(b_k) - f(a_k))(f(b_k) - f(c_k))} b_k + \frac{(y - f(a_k))(y - f(b_k))}{(f(c_k) - f(a_k))(f(c_k) - f(b_k))} c_k$$

Setting the target value $y = 0$ yields the unconstrained root estimate (Brent, 1973):

$$s = \frac{a_k f(b_k) f(c_k)}{(f(a_k) - f(b_k))(f(a_k) - f(c_k))} + \frac{b_k f(a_k) f(c_k)}{(f(b_k) - f(a_k))(f(b_k) - f(c_k))} + \frac{c_k f(a_k) f(b_k)}{(f(c_k) - f(a_k))(f(c_k) - f(b_k))}$$

When two ordinate values coincide, IQI degenerates via division by zero, prompting the algorithm to switch to a linear secant interpolation step between $a_k$ and $b_k$ (Brent, 1973; Dekker, 1969):

$$s = b_k - f(b_k) \frac{b_k - a_k}{f(b_k) - f(a_k)}$$

To prevent slow asymptotic convergence or numerical instability, Brent introduced five bounding checks that dictate whether the interpolated candidate $s$ is accepted or rejected in favor of an exact interval bisection step $m = \frac{a_k + b_k}{2}$ (Brent, 1973):

1. The candidate point $s$ must fall strictly within the interior interval bounded by $\frac{3a_k + b_k}{4}$ and $b_k$.

2. If the previous iteration executed a bisection step, the proposed step length $|s - b_k|$ must be strictly less than half the magnitude of the step taken two iterations prior: $|s - b_k| < \frac{1}{2}|b_k - b_{k-1}|$.

3. If the previous iteration executed an interpolation step, the proposed step length must be strictly less than half the distance between the two preceding historical iterates: $|s - b_k| < \frac{1}{2}|b_{k-1} - b_{k-2}|$.

4. If the previous iteration executed a bisection step, the magnitude of the historical interval must exceed the current convergence tolerance: $|b_k - b_{k-1}| > \delta$, where $\delta = 2 \varepsilon_{\text{mach}} |b_k| + t_{\text{abs}}$.

5. If the previous iteration executed an interpolation step, the previous historical interval must also satisfy this criterion: $|b_{k-1} - b_{k-2}| > \delta$.

These safeguards ensure that whenever an interpolation step stalls or exhibits erratic behavior, a bisection step is forced, guaranteeing that the enclosing interval shrinks by at least a factor of two over successive steps (Brent, 1973).

### Breakdown Modes, Singularities, and Pathologies

Newton-type algorithms can fail under several distinct mathematical conditions:

* Derivative Singularities: The Jacobian operator $F'(x_k)$ becomes singular, meaning $\det(F'(x_k)) = 0$, or its condition number diverges to infinity (Polyak, 2007). In this limit, the affine model loses uniqueness, and the inverse operator $[F'(x_k)]^{-1}$ is unbounded.

* Violation of the Semi-Local Criterion: When the Kantorovich parameter satisfies $h = \beta \eta L > \frac{1}{2}$, the discriminant of the majorizing polynomial $p(t)$ becomes strictly negative (Kantorovich, 1948). The contraction mapping property of the step fails, leading to periodic limit cycles, chaotic orbits, or divergence (Polyak, 2007).

* Basin Boundary Instability: In nonlinear systems with multiple isolated zeros, the basins of attraction $\mathcal{B}(x^*_i) = \{x_0 \mid \lim_{k\to\infty} x_k = x^*_i\}$ form fractal boundaries (Polyak, 2007). Initial guesses placed near these Julia sets exhibit extreme sensitivity to roundoff errors, producing large pseudo-random displacements.

* Root Multiplicity Degeneracy: When evaluating a root of multiplicity $m > 1$, the first derivative vanishes at the solution: $F'(x^*) = 0$. The asymptotic rate of convergence drops from quadratic to linear, with an asymptotic error constant $C = 1 - \frac{1}{m}$.

Brent's method can also encounter specific failure modes:

* Discontinuous Functions: If $f(x)$ contains an essential, step, or pole discontinuity across which its sign alternates, the algorithm converges directly to the singularity rather than a root (Brent, 1973).

* Even Multiplicities: Roots with an even integer multiplicity ($f(x^*) = 0$ with $f(x) \ge 0$ locally) violate the sign change condition $f(a)f(b) \le 0$, preventing algorithm initialization unless the global search detects the minimum (Brent, 1973).

### Recovery Protocols and Algorithmic Fallbacks

When the pure Newton iteration encounters a singular or ill-conditioned Jacobian, several stabilization techniques can be applied:

* Levenberg-Marquardt Regularization: The singular linear system is replaced with a damped least-squares problem parameterized by a dynamic Tikhonov factor $\lambda_k > 0$ (Polyak, 2007):

$$(F'(x_k)^T F'(x_k) + \lambda_k I) s_k = -F'(x_k)^T F(x_k)$$

As $\lambda_k \to \infty$, the search vector shifts toward the steepest descent direction $s_k \to -\frac{1}{\lambda_k} F'(x_k)^T F(x_k)$, restoring well-posedness when the Jacobian loses rank (Polyak, 2007).

* Armijo Line Search Backtracking: The full step is scaled by an adaptive parameter $\alpha_k \in (0, 1]$ satisfying the sufficient decrease condition on the merit function $\theta(x) = \frac{1}{2} \Vert{}F(x)\Vert{}_2^2$:

$$\theta(x_k + \alpha_k s_k) \le \theta(x_k) + c_1 \alpha_k \nabla \theta(x_k)^T s_k$$

* Schroeder Acceleration for Multiple Roots: If the multiplicity $m \ge 2$ is known, quadratic convergence is restored using the modified step $x_{k+1} = x_k - m [F'(x_k)]^{-1} F(x_k)$. When $m$ is unknown, it can be estimated dynamically from the ratio of successive steps:

$$m \approx \left(1 - \frac{\|x_{k+1} - x_k\|}{\|x_k - x_{k-1}\|}\right)^{-1}$$

In Brent's method, if a proposed interpolation step fails any of the five safeguard checks, the algorithm immediately executes a bisection step (Brent, 1973). If roundoff errors prevent the brackets from narrowing before the function values reach zero, the algorithm terminates cleanly using the floating-point tolerance check:

$$|b_k - a_k| \le 2 \varepsilon_{\text{mach}} |b_k| + t_{\text{abs}}$$

### Precision Boundaries and Numerical Stability Limits

In standard double-precision floating-point arithmetic (IEEE 754 binary64, where machine epsilon $\varepsilon_{\text{mach}} = 2^{-52} \approx 2.22 \times 10^{-16}$), computing the residual $F(x_k)$ is subject to roundoff error. The absolute limit of accuracy for a computed root $\hat{x}^*$ is dictated by the condition number of the root-finding problem:

$$\text{cond}(F, x^*) = \|[F'(x^*)]^{-1}\|_{\mathcal{L}(Y, X)}$$

The limiting forward error is bounded by:

$$\|x^* - \hat{x}^*\| \le \|[F'(x^*)]^{-1}\| \cdot \varepsilon_{\text{mach}} \|F(x^*)\| + \mathcal{O}(\varepsilon_{\text{mach}}^2)$$

If the Jacobian condition number satisfies $\kappa(F'(x_k)) = \|F'(x_k)\| \|[F'(x_k)]^{-1}\| \ge \varepsilon_{\text{mach}}^{-1}$, solving the linear system introduces catastrophic cancellation, which corrupts the update direction.

For Brent's method, the search interval cannot be contracted beyond the machine precision limit:

$$|b - a| \le \varepsilon_{\text{mach}} |b|$$

Below this threshold, the floating-point midpoint evaluates identically to one of the interval endpoints ($fl((a+b)/2) \in \{a, b\}$), causing the algorithm to stall in an infinite loop unless the step tolerance check is enforced (Brent, 1973).

## Unconstrained Nonlinear Optimization: The BFGS Quasi-Newton Framework

### Variational Foundations and Dual Space Formulations

Unconstrained optimization addresses the minimization of an objective function $f: \mathbb{R}^n \to \mathbb{R}$, where $f \in \mathcal{C}^2$ (Nocedal & Wright, 2006). Classical Newton optimization uses the inverse Hessian $[\nabla^2 f(x_k)]^{-1}$ to adjust the search direction based on local curvature:

$$p_k = -[\nabla^2 f(x_k)]^{-1} \nabla f(x_k)$$

Because assembling and factorizing the exact Hessian requires $\mathcal{O}(n^3)$ operations and does not guarantee positive definiteness in non-convex regions, quasi-Newton methods maintain a symmetric positive definite approximation $B_k \approx \nabla^2 f(x_k)$, or its direct inverse $H_k = B_k^{-1}$ (Broyden, 1970; Fletcher, 1970; Goldfarb, 1970; Shanno, 1970).

Let the spatial displacement vector be $s_k = x_{k+1} - x_k = \alpha_k p_k$, and let the gradient variation vector be $y_k = \nabla f(x_{k+1}) - \nabla f(x_k)$. The updated operator $B_{k+1}$ must satisfy the secant equation (Nocedal & Wright, 2006):

$$B_{k+1} s_k = y_k \quad \Longleftrightarrow \quad H_{k+1} y_k = s_k$$

The Broyden-Fletcher-Goldfarb-Shanno (BFGS) update is derived variationally by finding the minimal symmetric perturbation from $B_k$ that satisfies the secant condition, measured under a weighted Frobenius norm (Goldfarb, 1970; Nocedal & Wright, 2006):

$$\min_{B} \|W^{-1/2} (B - B_k) W^{-1/2}\|_F^2 \quad \text{subject to} \quad B = B^T, \quad B s_k = y_k$$

Choosing the weighting matrix as the average Hessian $W = \int_0^1 \nabla^2 f(x_k + \tau s_k) \, d\tau$ yields the rank-two update:

$$B_{k+1} = B_k - \frac{B_k s_k s_k^T B_k}{s_k^T B_k s_k} + \frac{y_k y_k^T}{y_k^T s_k}$$

Applying the Sherman-Morrison-Woodbury inversion theorem produces the inverse update formula, which computes $H_{k+1} \approx [\nabla^2 f(x_{k+1})]^{-1}$ in $\mathcal{O}(n^2)$ operations (Nocedal & Wright, 2006):

$$H_{k+1} = (I - \rho_k s_k y_k^T) H_k (I - \rho_k y_k s_k^T) + \rho_k s_k s_k^T, \quad \text{where} \quad \rho_k = \frac{1}{y_k^T s_k}$$

### Curvature Preservation and the Wolfe Invariant

To ensure that every quasi-Newton direction $p_{k+1} = -H_{k+1} \nabla f(x_{k+1})$ remains a valid descent direction ($\nabla f(x_{k+1})^T p_{k+1} < 0$), the operator $H_{k+1}$ must remain strictly positive definite (Nocedal & Wright, 2006). This property holds if and only if the curvature condition is satisfied:

$$y_k^T s_k > 0$$

Theorem (Preservation of Positive Definiteness):

Let $H_k$ be a real, symmetric, positive definite matrix ($H_k \succ 0$). If the spatial increment $s_k$ and the gradient change $y_k$ satisfy the curvature inequality $y_k^T s_k > 0$, then the updated matrix $H_{k+1}$ generated by the BFGS formula is symmetric and strictly positive definite ($H_{k+1} \succ 0$) (Nocedal & Wright, 2006).

Proof:

Symmetry is evident from the algebraic structure of the update equation. Let $z \in \mathbb{R}^n$ be any non-zero vector ($z \ne 0$). The associated quadratic form is:

$$z^T H_{k+1} z = z^T (I - \rho_k s_k y_k^T) H_k (I - \rho_k y_k s_k^T) z + \rho_k (z^T s_k)^2$$

Define the transformed vector $v = (I - \rho_k y_k s_k^T) z$. Substituting $v$ into the quadratic form gives:

$$z^T H_{k+1} z = v^T H_k v + \rho_k (z^T s_k)^2$$

Case 1: Suppose $v \ne 0$. Because $H_k$ is strictly positive definite by assumption, $v^T H_k v > 0$. Furthermore, because $y_k^T s_k > 0$, the scalar factor is positive: $\rho_k = (y_k^T s_k)^{-1} > 0$. The squared inner product term satisfies $\rho_k (z^T s_k)^2 \ge 0$, which ensures $z^T H_{k+1} z > 0$.

Case 2: Suppose $v = 0$. By definition:

$$(I - \rho_k y_k s_k^T) z = 0 \implies z = \rho_k (s_k^T z) y_k$$

Taking the inner product of both sides with $s_k$ yields:

$$s_k^T z = \rho_k (s_k^T z) (s_k^T y_k) = (s_k^T z) \frac{s_k^T y_k}{y_k^T s_k} = s_k^T z$$

This identity holds identically. However, if $s_k^T z = 0$, then $z = \rho_k (0) y_k = 0$, which contradicts the initial condition that $z \ne 0$. Therefore, $s_k^T z$ must be nonzero ($s_k^T z \ne 0$). Under these conditions, the quadratic form evaluates to:

$$z^T H_{k+1} z = 0 + \rho_k (s_k^T z)^2 > 0$$

Thus, $z^T H_{k+1} z > 0$ for all nonzero vectors $z \in \mathbb{R}^n \setminus \{0\}$, confirming that $H_{k+1} \succ 0$. $\blacksquare$

To enforce $y_k^T s_k > 0$ automatically during numerical line search, the step length $\alpha_k$ must satisfy the Wolfe conditions (Nocedal & Wright, 2006):

1. Armijo Sufficient Decrease:

$$f(x_k + \alpha_k p_k) \le f(x_k) + c_1 \alpha_k \nabla f(x_k)^T p_k, \quad 0 < c_1 < 1$$

2. Curvature Condition:

$$\nabla f(x_k + \alpha_k p_k)^T p_k \ge c_2 \nabla f(x_k)^T p_k, \quad c_1 < c_2 < 1$$

Subtracting $\nabla f(x_k)^T p_k$ from both sides of the curvature condition gives:

$$(\nabla f(x_k + \alpha_k p_k) - \nabla f(x_k))^T p_k \ge (c_2 - 1) \nabla f(x_k)^T p_k$$

Multiplying by $\alpha_k > 0$ and substituting $s_k = \alpha_k p_k$ yields:

$$y_k^T s_k \ge (c_2 - 1) \alpha_k \nabla f(x_k)^T p_k$$

Because $p_k$ is a descent direction, $\nabla f(x_k)^T p_k < 0$. Because $c_2 < 1$, the right-hand factor $(c_2 - 1)$ is strictly negative. The product of these two negative terms is strictly positive:

$$y_k^T s_k > 0$$

The Wolfe curvature condition thus enforces positive definiteness at every update step (Nocedal & Wright, 2006).

### Breakdown Modes, Indefiniteness, and Remediation Protocols

The BFGS algorithm can degrade or fail under several specific conditions:

* Non-Convexity and Negative Curvature: In non-convex regions, if an unconstrained or poorly configured line search fails to satisfy the Wolfe conditions, it can generate an update where $y_k^T s_k \le 0$, destroying positive definiteness.

* Roundoff Accumulation and Loss of Symmetry: In finite-precision floating-point arithmetic, the outer-product terms in the rank-two update accumulate asymmetric roundoff errors over many iterations, eventually causing $H_k$ to develop negative eigenvalues.

* Memory Saturation: For large optimization problems with dimension $n \ge 10^5$, storing and updating the dense $n \times n$ matrix $H_k$ requires $\mathcal{O}(n^2)$ memory, leading to memory exhaustion.

These issues are resolved through specific algorithmic modifications:

* Powell's Damped BFGS Update: When $y_k^T s_k$ becomes too small, $y_k$ is modified to preserve positive definiteness (Powell, 1978):

$$r_k = \theta_k y_k + (1 - \theta_k) B_k s_k$$

where the interpolation weight $\theta_k \in [0, 1]$ is defined by:

$$\theta_k = \begin{cases} 1 & \text{if } y_k^T s_k \ge 0.2 \, s_k^T B_k s_k \\ \frac{0.8 \, s_k^T B_k s_k}{s_k^T B_k s_k - y_k^T s_k} & \text{if } y_k^T s_k < 0.2 \, s_k^T B_k s_k \end{cases}$$

Replacing $y_k$ with $r_k$ in the update guarantees that $s_k^T r_k \ge 0.2 \, s_k^T B_k s_k > 0$.

* Explicit Symmetrization: Roundoff-induced asymmetry is eliminated by projecting $H_k$ back onto the symmetric subspace at the end of each iteration:

$$H_k \leftarrow \frac{1}{2} (H_k + H_k^T)$$

* Limited-Memory BFGS (L-BFGS): Rather than forming $H_k$ explicitly, the L-BFGS variant retains only the $m$ most recent vector pairs $\{s_i, y_i\}_{i=k-m}^{k-1}$ (typically $m \in [5, 30]$). The product $H_k \nabla f(x_k)$ is evaluated via the two-loop recursion algorithm in $\mathcal{O}(mn)$ operations and memory (Nocedal & Wright, 2006).

* Descent Safeguard and Directional Reset: If the computed search direction fails the descent condition:

$$p_k^T \nabla f(x_k) \ge -\varepsilon_{\text{mach}} \|p_k\|_2 \|\nabla f(x_k)\|_2$$

the current inverse Hessian approximation is discarded and reset to the identity matrix ($H_k \leftarrow I$), reverting to a steepest descent step.

### Condition Limits and Superlinear Stagnation

Under the standard Dennis-Moré characterization, the BFGS sequence converges superlinearly if and only if the search direction asymptotically matches the true Newton step (Nocedal & Wright, 2006):

$$\lim_{k \to \infty} \frac{\|(B_k - \nabla^2 f(x^*)) p_k\|}{\|p_k\|} = 0$$

However, the condition number of the Hessian approximation tends to grow across iterations:

$$\kappa(B_k) = \|B_k\|_2 \|B_k^{-1}\|_2 \approx \kappa(\nabla^2 f(x^*))$$

When $\kappa(B_k) \ge \varepsilon_{\text{mach}}^{-1} \approx 10^{16}$, the search direction vector $p_k = -H_k \nabla f(x_k)$ becomes numerically orthogonal to the true gradient, causing the line search to fail.

## Direct Linear Solvers: LU Decomposition and Gaussian Elimination

### Error Mechanics and the Gaussian Growth Factor

Gaussian Elimination with Partial Pivoting (GEPP) computes the factorization of a row-permuted square matrix $A \in \mathbb{R}^{n \times n}$ into a unit lower triangular matrix $L$ and an upper triangular matrix $U$ (Higham, 2002):

$$P A = L U$$

At each step $k$, row permutations ensure that the pivot element satisfies the partial pivoting condition:

$$|a_{kk}^{(k)}| = \max_{i \ge k} |a_{ik}^{(k)}|$$

The subdiagonal elimination multipliers $\ell_{ik} = \frac{a_{ik}^{(k)}}{a_{kk}^{(k)}}$ are bounded by construction:

$$|\ell_{ik}| \le 1, \quad \forall i > k$$

The active submatrix elements are updated via outer products:

$$a_{ij}^{(k+1)} = a_{ij}^{(k)} - \ell_{ik} a_{kj}^{(k)}, \quad \forall i, j \ge k+1$$

In finite-precision arithmetic, Wilkinson's backward error analysis demonstrates that the computed factors $\hat{L}$ and $\hat{U}$ satisfy the perturbed system (Higham, 2002; Wilkinson, 1961):

$$\hat{L} \hat{U} = P A + \Delta A$$

where the backward perturbation matrix $\Delta A$ satisfies the componentwise bound (Higham, 2002):

$$|\Delta A| \le \frac{n \varepsilon_{\text{mach}}}{1 - n \varepsilon_{\text{mach}}} |\hat{L}| |\hat{U}|$$

Taking matrix infinity-norms yields the standard backward error estimate (Higham, 2002):

$$\|\Delta A\|_\infty \le \gamma_n \rho_n(A) \|A\|_\infty$$

where $\gamma_n = \frac{n \varepsilon_{\text{mach}}}{1 - n \varepsilon_{\text{mach}}} = n \varepsilon_{\text{mach}} + \mathcal{O}(\varepsilon_{\text{mach}}^2)$, and $\rho_n(A)$ is the element growth factor (Higham, 2002; Wilkinson, 1961):

$$\rho_n(A) = \frac{\max_{i, j, k} |a_{ij}^{(k)}|}{\max_{i, j} |a_{ij}^{(1)}|}$$

Because $|\ell_{ik}| \le 1$, the magnitude of updated submatrix elements is bounded by:

$$|a_{ij}^{(k+1)}| \le |a_{ij}^{(k)}| + |\ell_{ik}| |a_{kj}^{(k)}| \le |a_{ij}^{(k)}| + |a_{kj}^{(k)}| \le 2 \max_{r, c} |a_{rc}^{(k)}|$$

Applying this bound inductively across $n - 1$ elimination stages yields the theoretical upper bound on element growth under partial pivoting (Wilkinson, 1961):

$$\rho_n(A) \le 2^{n-1}$$

### Pathological Instability: The Wilkinson Matrix

The theoretical upper bound $\rho_n(A) = 2^{n-1}$ is tight and is achieved by the Wilkinson matrix (Higham, Higham, & Pranesh, 2021; Wilkinson, 1961):

$$W_n = \begin{bmatrix}1 & 0 & 0 & \cdots & 0 & 1 \\-1 & 1 & 0 & \cdots & 0 & 1 \\-1 & -1 & 1 & \cdots & 0 & 1 \\\vdots & \vdots & \vdots & \ddots & \vdots & \vdots \\-1 & -1 & -1 & \cdots & 1 & 1 \\-1 & -1 & -1 & \cdots & -1 & 1\end{bmatrix}$$

Because the magnitude of the diagonal element in the active column is already tied for the maximum value ($\vert{}a_{kk}\vert{} = 1$), partial pivoting triggers no row interchanges ($P = I$) (Wilkinson, 1961). During elimination, the entries in the final column double at each stage:

$$a_{i, n}^{(k)} = 2^{k-1}, \quad \forall i \ge k$$

At the final elimination step $k = n$, the entry reaches $a_{nn}^{(n)} = 2^{n-1}$ (Wilkinson, 1961).

For dimension $n = 64$, the growth factor evaluates to $\rho_{64} = 2^{63} \approx 9.22 \times 10^{18}$. In IEEE double precision ($\varepsilon_{\text{mach}} \approx 2.22 \times 10^{-16}$), the backward error satisfies:

$$\|\Delta A\|_\infty \approx n \varepsilon_{\text{mach}} 2^{n-1} \|A\|_\infty \approx 64 \cdot (2.22 \times 10^{-16}) \cdot (9.22 \times 10^{18}) \|A\|_\infty \approx 1.31 \times 10^5 \|A\|_\infty$$

The backward perturbation exceeds the original matrix norm by five orders of magnitude. The computed solution $\hat{x}$ loses all significant digits to catastrophic cancellation, even though $W_n$ is well-conditioned ($\kappa_\infty(W_n) \approx n 2^{n-1}$ is modest for moderate $n$).

### Pivoting Paradigms, Condition Boundaries, and Remediation

When partial pivoting produces excessive element growth, alternative pivoting architectures must be deployed:

* Complete Pivoting (GECP): Permutes both rows and columns to select the globally maximal entry in the active submatrix as the pivot:

$$|a_{p, q}^{(k)}| = \max_{i, j \ge k} |a_{ij}^{(k)}| \implies P A Q = L U$$

Wilkinson proved that complete pivoting significantly restricts element growth (Wilkinson, 1961):

$$\rho_n^{CP} \le \sqrt{n} \left( 2^1 3^{1/2} 4^{1/3} \cdots n^{1/(n-1)} \right)^{1/2} \sim \mathcal{O}(n^{\frac{1}{4} \log n})$$

This stability comes at the expense of an increased search cost of $\mathcal{O}(n^3)$ comparisons, compared to $\mathcal{O}(n^2)$ for partial pivoting.

* Rook Pivoting: An intermediate strategy where the pivot element is chosen to be simultaneously the maximum entry in its column and its row:

$$|a_{rs}^{(k)}| \ge \max_{i \ge k} |a_{is}^{(k)}| \quad \text{and} \quad |a_{rs}^{(k)}| \ge \max_{j \ge k} |a_{rj}^{(k)}|$$

Rook pivoting limits element growth to small polynomial bounds while running in $\mathcal{O}(n^2)$ average comparisons (Higham, 2002).

* Householder QR Factorization: For problems where element growth cannot be permitted, the system $A x = b$ is solved using orthogonal transformations:

$$Q A = R \implies R x = Q^T b$$

Because orthogonal transformations preserve the Euclidean norm ($\Vert{}Q\Vert{}_2 = 1$), the growth factor is strictly bounded (Higham, 2002):

$$\rho_n^{QR} = 1$$

Householder QR factorization is unconditionally backward stable, satisfying $\|\Delta A\|_F \le c n \varepsilon_{\text{mach}} \|A\|_F$.

The forward error of the computed solution $\hat{x}$ across all direct linear solvers is governed by the condition number $\kappa(A) = \|A\| \|A^{-1}\|$:

$$\frac{\|x - \hat{x}\|}{\|x\|} \le \frac{\kappa(A) \frac{\|\Delta A\|}{\|A\|}}{1 - \kappa(A) \frac{\|\Delta A\|}{\|A\|}} \approx \kappa(A) \cdot n \cdot \rho_n(A) \cdot \varepsilon_{\text{mach}}$$

Whenever $\kappa(A) \cdot \rho_n(A) \ge \varepsilon_{\text{mach}}^{-1}$, all precision is lost.

In this regime, accuracy can be restored through mixed-precision iterative refinement. The residual vector is computed in extended precision:

$$r_k = b - A \hat{x}_k$$

The correction equation $A d_k = r_k$ is solved using the existing LU factors, and the solution is updated:

$$\hat{x}_{k+1} = \hat{x}_k + d_k$$

Iterative refinement restores forward accuracy to the machine precision level, provided $\kappa(A) \varepsilon_{\text{mach}} < 1$ (Higham, 2002).

## Large-Scale Spectral and Krylov Subspace Computations: The Lanczos Iteration

### Tridiagonalization Mechanics and Kaniel-Paige Convergence Theory

The Lanczos algorithm computes the extremal eigenvalues and eigenvectors of a large, sparse, real symmetric matrix $A = A^T \in \mathbb{R}^{n \times n}$ (Lanczos, 1950; Parlett, 1998). Initialized with a unit vector $v_1$ ($\Vert{}v_1\Vert{}_2 = 1$), the method constructs an orthonormal basis for the Krylov subspace:

$$\mathcal{K}_m(A, v_1) = \text{span}\{v_1, A v_1, A^2 v_1, \dots, A^{m-1} v_1\}$$

Because $A$ is symmetric, the Gram-Schmidt orthogonalization process reduces to a three-term recurrence (Lanczos, 1950; Parlett, 1998):

$$\beta_{j+1} v_{j+1} = w_j = A v_j - \alpha_j v_j - \beta_j v_{j-1}$$

where $v_0 = 0$, and the scalar coefficients are determined by the orthogonality conditions $v_j^T v_{j+1} = 0$ and $v_{j-1}^T v_{j+1} = 0$:

$$\alpha_j = v_j^T A v_j, \quad \beta_{j+1} = \|w_j\|_2$$

In matrix form, after $m$ iterations, the recurrence reads (Parlett, 1998):

$$A V_m = V_m T_m + \beta_{m+1} v_{m+1} e_m^T$$

where $V_m = [v_1, \dots, v_m] \in \mathbb{R}^{n \times m}$ has orthonormal columns ($V_m^T V_m = I_m$), $e_m = [0, \dots, 0, 1]^T$, and $T_m$ is a symmetric tridiagonal matrix:

$$T_m = \begin{bmatrix}\alpha_1 & \beta_2 & & \\\beta_2 & \alpha_2 & \ddots & \\& \ddots & \ddots & \beta_m \\& & \beta_m & \alpha_m\end{bmatrix}$$

The eigenvalues $\theta_1^{(m)} > \theta_2^{(m)} > \dots > \theta_m^{(m)}$ of $T_m$ are the Ritz values, and the vectors $x_i = V_m y_i$ (where $T_m y_i = \theta_i y_i$ and $\|y_i\|_2 = 1$) are the Ritz vectors (Parlett, 1998).

The convergence of Ritz values to the true eigenvalues $\lambda_1 \ge \lambda_2 \ge \dots \ge \lambda_n$ of $A$ is governed by the Kaniel-Paige-Saad convergence theory, which expresses the projection error through Chebyshev polynomials (Kaniel, 1966; Paige, 1971; Saad, 1980).

Theorem (Kaniel-Paige Bounds):

Let $\lambda_1$ be the maximal eigenvalue of $A$, and let $\phi_1 = \angle(v_1, u_1)$ denote the angle between the initial Lanczos vector $v_1$ and the normalized eigenvector $u_1$ corresponding to $\lambda_1$ (Kaniel, 1966; Saad, 1980). The error in the maximal Ritz value $\theta_1^{(m)}$ after $m$ iterations satisfies:

$$0 \le \lambda_1 - \theta_1^{(m)} \le (\lambda_1 - \lambda_n) \left( \frac{\tan \phi_1}{C_{m-1}(1 + 2 \gamma_1)} \right)^2$$

where $C_{m-1}(x)$ is the Chebyshev polynomial of the first kind of degree $m-1$, and $\gamma_1$ is the normalized spectral gap (Kaniel, 1966):

$$\gamma_1 = \frac{\lambda_1 - \lambda_2}{\lambda_2 - \lambda_n}$$

Using the explicit representation $C_{m-1}(1 + 2\gamma_1) = \frac{1}{2} \left[ R_1^{m-1} + R_1^{-(m-1)} \right]$, where $R_1 = 1 + 2\gamma_1 + 2\sqrt{\gamma_1^2 + \gamma_1} > 1$, the convergence bound simplifies to (Saad, 1980):

$$\lambda_1 - \theta_1^{(m)} \le 4 (\lambda_1 - \lambda_n) \tan^2(\phi_1) \cdot R_1^{-2(m-1)}$$

The Ritz values converge exponentially to isolated extreme eigenvalues, with the rate determined by the spectral gap ratio $\gamma_1$ (Kaniel, 1966; Saad, 1980).

### Finite Precision and Paige's Loss of Orthogonality Theorem

In exact arithmetic, the Lanczos basis vectors remain mutually orthogonal: $V_m^T V_m = I_m$. In finite-precision floating-point arithmetic, orthogonality degrades rapidly (Paige, 1971, 1976). The mechanics of this degradation were resolved by Christopher Paige (1971, 1976, 1980), who established that the loss of orthogonality is an analytical consequence of convergence.

Let computed quantities be denoted with hats. The finite-precision Lanczos recurrence satisfies (Paige, 1976):

$$A \hat{V}_m = \hat{V}_m \hat{T}_m + \hat{\beta}_{m+1} \hat{v}_{m+1} e_m^T + E_m, \quad \|E_m\|_2 \le c n \varepsilon_{\text{mach}} \|A\|_2$$

Paige proved that the loss of orthogonality does not arise from the uniform accumulation of roundoff errors across the basis (Paige, 1976, 1980). Instead, the basis vectors lose orthogonality along the direction of a converged Ritz vector.

Theorem (Paige's Loss of Orthogonality Bound):

Let $(\theta_j^{(m)}, y_j^{(m)})$ be an eigenpair of the computed tridiagonal matrix $\hat{T}_m$, and let $\hat{x}_j = \hat{V}_m y_j^{(m)}$ be the corresponding computed Ritz vector (Paige, 1976). The projection of the next Lanczos vector $\hat{v}_{m+1}$ onto the Ritz vector $\hat{x}_j$ satisfies:

$$|\hat{v}_{m+1}^T \hat{x}_j| = \frac{\tau_{m, j} \varepsilon_{\text{mach}} \|A\|_2}{\hat{\beta}_{m+1} |e_m^T y_j^{(m)}|}$$

where $\tau_{m, j} = \mathcal{O}(m)$ is a modest constant (Paige, 1976, 1980).

The denominator represents the residual norm of the computed Ritz pair:

$$\text{Residual}_j = \|A \hat{x}_j - \theta_j^{(m)} \hat{x}_j\|_2 = \hat{\beta}_{m+1} |e_m^T y_j^{(m)}|$$

As the Ritz value converges to a true eigenvalue of $A$, its residual norm approaches machine precision ($\text{Residual}_j \to \varepsilon_{\text{mach}} \Vert{}A\Vert{}_2$) (Paige, 1976). By Paige's theorem, this causes the projection $|\hat{v}_{m+1}^T \hat{x}_j|$ to approach unity:

$$|\hat{v}_{m+1}^T \hat{x}_j| \approx \frac{\varepsilon_{\text{mach}} \|A\|_2}{\varepsilon_{\text{mach}} \|A\|_2} = 1$$

The newly generated vector $\hat{v}_{m+1}$ loses linear independence and begins aligning with the already converged Ritz vector $\hat{x}_j$ (Paige, 1976, 1980).

### Ghost Eigenvalues and Reorthogonalization Architectures

This alignment causes the algorithm to restart a secondary, spurious Krylov recurrence along that eigenvector (Paige, 1976; Parlett, 1998). As a result, duplicate copies of the converged eigenvalue—known as ghost eigenvalues—appear in the spectrum of the tridiagonal matrix $\hat{T}_m$ (Cullum & Willoughby, 1985; Parlett & Scott, 1979). These ghost eigenvalues do not reflect true multiplicity in $A$; they are numerical artifacts of lost orthogonality.

Four main architectures are used to stabilize the Lanczos iteration:

* Full Reorthogonalization (FRO): At each step $j$, the newly generated vector $w_j$ is explicitly orthogonalized against all previous vectors $v_1, \dots, v_j$ using two passes of Modified Gram-Schmidt:

$$w_j \leftarrow w_j - (v_i^T w_j) v_i, \quad \forall i = 1, \dots, j$$

This maintains machine-level orthogonality ($\Vert{}V_m^T V_m - I\Vert{}_2 \le \mathcal{O}(m \varepsilon_{\text{mach}})$) but requires storing all basis vectors and incurs an $\mathcal{O}(m^2 n)$ computational cost.

* Selective Orthogonalization (SO - Parlett & Scott): Instead of orthogonalizing against all historical vectors, the algorithm tracks the Ritz residuals $\beta_{j+1} \vert{}e_j^T y_i^{(j)}\vert{}$. When a Ritz pair converges to within $\sqrt{\varepsilon_{\text{mach}}} \Vert{}A\Vert{}_2$, the next vector $v_{j+1}$ is orthogonalized exclusively against that specific Ritz vector $x_i = V_j y_i$. This maintains semi-orthogonality ($\Vert{}V_m^T V_m - I\Vert{}_2 \le \sqrt{\varepsilon_{\text{mach}}}$) with minimal overhead (Parlett & Scott, 1979).

* Partial Reorthogonalization (PRO - Simon): Maintains a scalar recurrence relation that bounds the inner products $\omega_{i, j} = v_i^T v_j$ in floating-point arithmetic. When $\max_i \vert{}\omega_{i, j+1}\vert{}$ exceeds a threshold $\tau \approx \sqrt{\varepsilon_{\text{mach}}}$, a full Gram-Schmidt reorthogonalization is triggered for that step (Simon, 1984).

* Cullum-Willoughby Identification Test: Executes the three-term recurrence without any reorthogonalization. The spectrum of $T_m$ is compared against the spectrum of the submatrix $\tilde{T}_{m-1}$ formed by deleting the first row and column of $T_m$. By the properties of Sturm sequences, any eigenvalue that appears in both $T_m$ and $\tilde{T}_{m-1}$ to within machine precision is identified as a ghost eigenvalue and discarded, isolating the true physical spectrum without storing the basis vectors (Cullum & Willoughby, 1985).

## Dense Eigenvalue and Singular Value Decompositions

### The Shifted QR Algorithm and the Wilkinson Shift

For dense, non-symmetric matrices, computing the Schur decomposition $A = Q T Q^T$ is carried out using the shifted QR algorithm (Francis, 1961; Wilkinson, 1965). The matrix is first reduced to upper Hessenberg form (or symmetric tridiagonal form if $A = A^T$) using Householder reflections (Wilkinson, 1965). The unshifted QR iteration factors $A_k = Q_k R_k$ and computes $A_{k+1} = R_k Q_k = Q_k^T A_k Q_k$, which preserves the Hessenberg structure (Francis, 1961).

To accelerate convergence, shifts $\mu_k$ are introduced:

$$\begin{aligned} A_k - \mu_k I &= Q_k R_k \\ A_{k+1} &= R_k Q_k + \mu_k I = Q_k^T A_k Q_k \end{aligned}$$

For symmetric tridiagonal matrices:

$$T = \begin{bmatrix} a_1 & b_1 & & \\ b_1 & a_2 & \ddots & \\ & \ddots & \ddots & b_{n-1} \\ & & b_{n-1} & a_n \end{bmatrix}$$

the Wilkinson shift selects the eigenvalue of the trailing $2 \times 2$ principal submatrix:

$$B_n = \begin{bmatrix} a_{n-1} & b_{n-1} \\ b_{n-1} & a_n \end{bmatrix}$$

that is closest to the corner element $a_n$ (Wilkinson, 1965). Let $d = \frac{a_{n-1} - a_n}{2}$. The shift is given by:

$$\mu = a_n - \frac{\text{sign}(d) b_{n-1}^2}{\vert{}d\vert{} + \sqrt{d^2 + b_{n-1}^2}}$$

where $\text{sign}(0) = 1$.

Convergence Rate and Asymptotic Behavior:

The Wilkinson shift guarantees global convergence for symmetric tridiagonal matrices and achieves asymptotic cubic convergence ($q = 3$) (Parlett, 1998; Wilkinson, 1965):

$$\vert{}b_{n-1}^{(k+1)}\vert{} \le \frac{\vert{}b_{n-1}^{(k)}\vert{}^3 \vert{}b_{n-2}^{(k)}\vert{}^2}{\vert{}d_k\vert{}^4} \implies \vert{}b_{n-1}^{(k+1)}\vert{} = \mathcal{O}\left(\vert{}b_{n-1}^{(k)}\vert{}^3\right)$$

The off-diagonal element $b_{n-1}$ vanishes within 2 to 3 iterations per eigenvalue, allowing deflation to isolate $a_n$ as an eigenvalue and reduce the active matrix dimension (Parlett, 1998).

For non-symmetric matrices with complex conjugate eigenvalues, a single real shift can stall convergence. The Francis Implicit Double-Shift algorithm avoids complex arithmetic by combining two conjugate shifts $\mu_1, \bar{\mu}_2$ into a single real degree-two polynomial step (Francis, 1961):

$$(A - \mu_1 I)(A - \bar{\mu}_2 I) = A^2 - 2 \text{Re}(\mu_1) A + \vert{}\mu_1\vert{}^2 I = Q R$$

By the Implicit Q Theorem, computing only the first column of this matrix product and chasing the resulting bulge down the Hessenberg form using orthogonal similarity transforms maintains real arithmetic while preserving quadratic convergence (Francis, 1961; Golik, 2004).

### Golub-Kahan Bidiagonalization and Demmel-Kahan High-Accuracy SVD

To compute the singular value decomposition (SVD) $A = U \Sigma V^T$ of a general matrix $A \in \mathbb{R}^{m \times n}$ ($m \ge n$), forming the product $A^T A$ explicitly degrades the condition number:

$$\kappa(A^T A) = (\kappa(A))^2$$

If $\kappa(A) \ge \varepsilon_{\text{mach}}^{-1/2} \approx 10^8$, singular values smaller than $\sqrt{\varepsilon_{\text{mach}}} \sigma_{\max}$ are lost to subtractive cancellation (Demmel & Kahan, 1990).

The Golub-Kahan algorithm avoids this by applying alternating Householder reflections from the left and right directly to $A$ (Golub & Kahan, 1965):

$$U_0^T A V_0 = B = \begin{bmatrix} \alpha_1 & \beta_1 & & \\ & \alpha_2 & \ddots & \\ & & \ddots & \beta_{n-1} \\ & & & \alpha_n \\ \hline & & 0 & \end{bmatrix}$$

Computing the singular values of $B$ is mathematically equivalent to computing the eigenvalues of $B^T B$ (Golub & Kahan, 1965). However, applying an explicit shift $\mu^2$ to $B^T B$ requires evaluating differences of squares, which introduces subtractive cancellation for small singular values (Demmel & Kahan, 1990).

The Demmel-Kahan Zero-Shift QR algorithm addresses this by setting the shift to zero ($\mu = 0$) (Demmel & Kahan, 1990). In this case, the QR transformation corresponds to an implicit Cholesky factorization:

$$B^T B = R^T Q^T Q R = R^T R \implies R = \text{Cholesky}(B^T B)$$

Demmel and Kahan showed that computing this step using Givens rotations from top to bottom requires no subtractions:

$$\begin{aligned} & \text{Initialize } c_0 = 1, \quad s_0 = 0 \\ & \text{For } i = 1 \text{ to } n-1: \\ & \quad \begin{bmatrix} c_i & s_i \\ -s_i & c_i \end{bmatrix} \begin{bmatrix} \alpha_i c_{i-1} \\ \beta_i \end{bmatrix} = \begin{bmatrix} \bar{\alpha}_i \\ 0 \end{bmatrix} \implies c_i = \frac{\alpha_i c_{i-1}}{\sqrt{(\alpha_i c_{i-1})^2 + \beta_i^2}}, \quad s_i = \frac{\beta_i}{\sqrt{(\alpha_i c_{i-1})^2 + \beta_i^2}} \\ & \quad \bar{\beta}_i = s_i \alpha_{i+1}, \quad \alpha_{i+1} \leftarrow c_i \alpha_{i+1} \end{aligned}$$

Because this update involves only products, quotients, and square roots of positive sums, it introduces no subtractive cancellation (Demmel & Kahan, 1990).

Theorem (Demmel-Kahan High Relative Accuracy):

Let $\sigma_1 \ge \sigma_2 \ge \dots \ge \sigma_n$ be the exact singular values of the bidiagonal matrix $B$, and let $\hat{\sigma}_i$ be the singular values computed by the zero-shift QR algorithm in floating-point arithmetic with machine precision $\varepsilon_{\text{mach}}$ (Demmel & Kahan, 1990). Every computed singular value satisfies the componentwise relative error bound:

$$\frac{\vert{}\hat{\sigma}_i - \sigma_i\vert{}}{\sigma_i} \le c \cdot n \cdot \varepsilon_{\text{mach}} + \mathcal{O}(\varepsilon_{\text{mach}}^2)$$

where $c$ is a small integer constant independent of the matrix condition number $\kappa(B)$.

This bound ensures that tiny singular values near the floating-point underflow threshold are computed with the same relative accuracy as the largest singular values (Demmel & Kahan, 1990).

### Deflation Dynamics and Eigenvector Sensitivity

The QR algorithm for eigenvalues and the SVD rely on deflation to decouple submatrices (Francis, 1961; Parlett, 1998). A subdiagonal element $b_j$ is set to zero when it satisfies the threshold criterion:

$$\vert{}b_j\vert{} \le \text{tol} \cdot (\vert{}a_j\vert{} + \vert{}a_{j+1}\vert{})$$

If the deflation threshold is set below machine precision ($\text{tol} < \varepsilon_{\text{mach}}$), roundoff prevents the condition from triggering, causing stagnation. If $\text{tol}$ is set too large, deflation introduces an uncontrolled perturbation of order $\vert{}b_j\vert{}$ into the spectrum.

For non-normal matrices, eigenvalue sensitivity is governed by the Bauer-Fike Theorem (Bauer & Fike, 1960).

Theorem (Bauer-Fike):

Let $A \in \mathbb{C}^{n \times n}$ be a diagonalizable matrix with spectral decomposition $A = X \Lambda X^{-1}$. If $\hat{\lambda}$ is an eigenvalue of a perturbed matrix $A + E$, then:

$$\min_{\lambda \in \Lambda(A)} \vert{}\hat{\lambda} - \lambda\vert{} \le \kappa_p(X) \Vert{}E\Vert{}_p$$

where $\kappa_p(X) = \Vert{}X\Vert{}_p \Vert{}X^{-1}\Vert{}_p$ is the condition number of the eigenvector matrix $X$.

When the eigenvectors are nearly linearly dependent ($\kappa_p(X) \gg 1$), small backward errors $\Vert{}E\Vert{} \sim \varepsilon_{\text{mach}}$ can produce large perturbations in the computed eigenvalues, even though the QR algorithm itself is backward stable (Bauer & Fike, 1960).

| **Algorithmic Class**          | **Representative Algorithm** | **Asymptotic Convergence Rate**                             | **Dominant Breakdown Mode**                                                   | **Fallback / Remediation Architecture**                                                  | **Reliability Bound Threshold**                 |
| ------------------------------ | ---------------------------- | ----------------------------------------------------------- | ----------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------- | ----------------------------------------------- |
| Nonlinear Root-Finding         | Newton-Kantorovich           | Quadratic (q = 2) (Argyros, 1998; Kantorovich, 1948)        | Derivative singularity or h > 0.5 (Kantorovich, 1948; Polyak, 2007)           | Levenberg-Marquardt regularization, Line search (Polyak, 2007)                           | cond(F'(x_k)) >= 1 / eps_mach                   |
| Bracketed Root-Finding         | Brent's Method               | Superlinear (q approx 1.839) (Brent, 1973)                  | Essential discontinuities or even multiplicities (Brent, 1973)                | Automatic fallback to bisection (Brent, 1973)                                            |                                                 |
| Unconstrained Optimization     | BFGS Quasi-Newton            | Superlinear (q > 1) (Nocedal & Wright, 2006)                | Curvature violation (y^T s <= 0) or loss of symmetry (Nocedal & Wright, 2006) | Powell damping, L-BFGS two-loop recursion (Nocedal & Wright, 2006; Powell, 1978)         | cond(B_k) >= 1 / eps_mach                       |
| Direct Linear Solvers          | GEPP (PA = LU)               | Finite ((2/3) n^3 operations) (Higham, 2002)                | Exponential element growth (rho_n -> 2^(n-1)) (Wilkinson, 1961)               | Rook pivoting, Complete pivoting, Householder QR (Higham, 2002; Wilkinson, 1961)         | cond(A) * rho_n(A) >= 1 / eps_mach              |
| Sparse Hermitian Eigensolvers  | Symmetric Lanczos            | Exponential (Chebyshev-bounded) (Kaniel, 1966; Paige, 1971) | Loss of orthogonality, Ghost eigenvalues (Paige, 1976, 1980)                  | Selective (SO) or Partial (PRO) reorthogonalization (Parlett & Scott, 1979; Simon, 1984) | Residual <= eps_mach *                          |
| Dense Tridiagonal Eigensolvers | QR with Wilkinson Shift      | Cubic (q = 3) (Parlett, 1998; Wilkinson, 1965)              | Stagnation on complex pairs (non-Hermitian) (Francis, 1961; Wilkinson, 1965)  | Francis Implicit Double-Shift, Deflation safeguards (Francis, 1961; Parlett, 1998)       | Subdiagonal                                     |
| Bidiagonal SVD Solvers         | Demmel-Kahan QR              | Linear (zero-shift Cholesky) (Demmel & Kahan, 1990)         | Slow convergence for large singular values (Golub & Kahan, 1965)              | Hybrid differential quotient-difference with shifts (dqds) (Fernando & Parlett, 1994)    | sigma_min <= Underflow Limit (approx 10^(-308)) |

## Synthesis: Error Regimes, Precision Limits, and Fallback Hierarchies

### Error Propagation and Perturbation Regimes

The behavior of numerical algorithms across finite-precision environments is governed by the separation between the condition number of the problem and the backward stability of the algorithm:

$$\text{Forward Error} \le \text{Condition Number} \times \text{Backward Error}$$

Numerical methods fall into three broad stability categories:

* Backward Stable Algorithms: Methods where the computed solution satisfies an exactly perturbed problem $(A + \Delta A)\hat{x} = b$ with $\Vert{}\Delta A\Vert{} / \Vert{}A\Vert{} \le c_n \varepsilon_{\text{mach}}$. Householder QR factorization, Givens rotations, and the shifted QR algorithm are backward stable (Higham, 2002; Wilkinson, 1965). For these methods, errors are driven entirely by the problem's conditioning, not by the accumulation of roundoff errors.

* Conditionally Stable Algorithms: Methods where backward stability depends on problem-specific parameters. In Gaussian elimination with partial pivoting, backward stability requires that the element growth factor satisfy $\rho_n(A) = \mathcal{O}(1)$ or $\mathcal{O}(n^{1/2})$ (Higham, 2002; Wilkinson, 1961). If $\rho_n(A)$ approaches its theoretical limit of $2^{n-1}$, backward stability is lost.

* Structurally Unstable Algorithms: Methods where finite-precision roundoff alters the underlying mathematical structure. In the three-term symmetric Lanczos algorithm, roundoff errors break the theoretical orthogonality of the basis vectors (Paige, 1976, 1980). Without active stabilization (such as reorthogonalization or ghost-filtering), the method cannot be run reliably over many iterations.

### Systematic Stabilization Hierarchies and Recovery Protocols

When an algorithm encounters numerical instability, recovery strategies should be deployed in a structured hierarchy based on the underlying failure mode:

When condition numbers approach the precision limit ($\kappa(A) \ge \varepsilon_{\text{mach}}^{-1}$), traditional solvers suffer from catastrophic cancellation. The primary remediation strategy is row and column equilibration. Applying diagonal scaling matrices $D_R A D_C \hat{x} = D_R b$ using the Ruiz or Sinkhorn-Knopp algorithms balances matrix norms and reduces the condition number. If the system remains ill-conditioned, computation should transition to mixed-precision iterative refinement (GMRES-IR) or higher-precision representations such as IEEE 754 binary128. For rank-deficient systems where $\sigma_{\min} \le \varepsilon_{\text{mach}} \sigma_{\max}$, standard inversion must be replaced by Tikhonov regularization:

$$x_\lambda = (A^T A + \lambda^2 I)^{-1} A^T b$$

where the parameter $\lambda$ is selected via the L-curve criterion to bound the effective condition number to $\kappa_{\text{eff}} \le \frac{\sigma_{\max}}{\lambda}$.

When Gaussian elimination exhibits large element growth ($\rho_n \gg 1$), partial pivoting should be aborted in favor of rook pivoting or complete pivoting (Higham, 2002; Wilkinson, 1961). If the active submatrix entries continue to grow, the linear solve should fall back to Householder QR factorization, which guarantees a growth factor of $\rho_n = 1$.

When the Lanczos algorithm experiences loss of orthogonality, tracking the Ritz residuals $\beta_{j+1} \vert{}e_j^T y_i\vert{}$ identifies when a Ritz pair converges to within $\sqrt{\varepsilon_{\text{mach}}} \Vert{}A\Vert{}_2$ (Paige, 1976; Parlett & Scott, 1979). At this threshold, selective reorthogonalization (SO) against the converged Ritz vectors prevents basis corruption. If memory limits prevent storing previous vectors, the Cullum-Willoughby test can be used to identify and remove ghost eigenvalues by comparing the spectrum of $T_m$ to that of its submatrix $\tilde{T}_{m-1}$ (Cullum & Willoughby, 1985).

When quasi-Newton optimization encounters non-positive curvature ($y_k^T s_k \le 0$), Powell's damped update should be engaged to scale the gradient change vector $y_k$ toward $B_k s_k$ (Powell, 1978). If the search direction fails to satisfy the descent property, the Hessian approximation should be reset to the identity matrix ($H_k \leftarrow I$), restarting the optimization along the steepest descent direction with an Armijo line search.

Finally, for problems with closely spaced eigenvalues or singular values that cannot be resolved within standard floating-point precision ($\Delta \lambda < \varepsilon_{\text{mach}} \Vert{}A\Vert{}$), computations should transition to arbitrary-precision arithmetic engines (such as MPFR) or interval arithmetic. By maintaining certified upper and lower bounds $[x_{\text{lower}}, x_{\text{upper}}]$ using directed rounding modes, interval methods provide guaranteed enclosures for the true mathematical solutions.

## References

Argyros, I. K. (1998). A new Kantorovich-type theorem for Newton's method. _Zeszyty Naukowe Politechniki Rzeszowskiej. Matematyka_, 26(23), 151–159.

Bauer, F. L., & Fike, C. T. (1960). Norms and exclusion theorems. _Numerische Mathematik_, 2(1), 137–141.

Brent, R. P. (1973). _Algorithms for Minimization without Derivatives_. Prentice-Hall.

Broyden, C. G. (1970). The convergence of a class of double-rank minimization algorithms 1. General considerations. _IMA Journal of Applied Mathematics_, 6(1), 76–90.

Ciarlet, P. G., & Mardare, C. (2014). On the Newton-Kantorovich theorem. _Mathematical Modelling and Numerical Analysis_, 48(4), 1181–1193.

Cullum, J. K., & Willoughby, R. A. (1985). _Lanczos Algorithms for Large Symmetric Eigenvalue Computations: Vol. 1: Theory_. Birkhäuser.

Dekker, T. J. (1969). Finding a zero by means of successive linear interpolation. In B. Dejon & P. Henrici (Eds.), _Constructive Aspects of the Fundamental Theorem of Algebra_ (pp. 37–48). Wiley-Interscience.

Demmel, J., & Kahan, W. (1990). Accurate singular values of bidiagonal matrices. _SIAM Journal on Scientific and Statistical Computing_, 11(5), 873–912.

Fernando, K. V., & Parlett, B. N. (1994). Accurate singular values and differential qd algorithms. _Numerische Mathematik_, 67(2), 191–229.

Fletcher, R. (1970). A new approach to variable metric algorithms. _The Computer Journal_, 13(3), 317–322.

Francis, J. G. F. (1961). The QR transformation a unitary analogue to the LR transformation—Part 1. _The Computer Journal_, 4(3), 265–271.

Goldfarb, D. (1970). A family of variable-metric methods derived by variational means. _Mathematics of Computation_, 24(109), 23–26.

Golik, W. L. (2004). The QR algorithm. _Eigenvalues of Dense Matrices_, Lecture Notes, 1–18.

Golub, G. H., & Kahan, W. (1965). Calculating the singular values and pseudo-inverse of a matrix. _Journal of the Society for Industrial and Applied Mathematics, Series B: Numerical Analysis_, 2(2), 205–224.

Higham, D. J., Higham, N. J., & Pranesh, S. (2021). Random matrices generating large growth in LU factorization with pivoting. _SIAM Journal on Matrix Analysis and Applications_, 42(1), 185–201.

Higham, N. J. (2002). _Accuracy and Stability of Numerical Algorithms_ (2nd ed.). Society for Industrial and Applied Mathematics.

Kaniel, S. (1966). Estimates for some executions of approximation of eigenvalues of symmetric operators. _Mathematics of Computation_, 20(93), 92–101.

Kantorovich, L. V. (1948). On Functional Equations and Functional Analysis. _Doklady Akademii Nauk SSSR_, 59(2), 209–212.

Kantorovich, L. V., & Akilov, G. P. (1982). _Functional Analysis_ (2nd ed.). Pergamon Press.

Lanczos, C. (1950). An iteration method for the solution of the eigenvalue problem of linear differential and integral operators. _Journal of Research of the National Bureau of Standards_, 45(4), 255–282.

Magreñán, A. (2016). On the convergence of Newton's method. _CMMSE Conference Proceedings_, 114–118.

Nocedal, J., & Wright, S. J. (2006). _Numerical Optimization_ (2nd ed.). Springer.

Ortega, J. M., & Rheinboldt, W. C. (1970). _Iterative Solution of Nonlinear Equations in Several Variables_. Academic Press.

Paige, C. C. (1971). _The Computation of Eigenvalues and Eigenvectors of Very Large Sparse Matrices_ (Doctoral dissertation). University of London.

Paige, C. C. (1976). Error analysis of the Lanczos algorithm for symmetric matrices. _IMA Journal of Applied Mathematics_, 18(3), 341–349.

Paige, C. C. (1980). Accuracy and effectiveness of the Lanczos algorithm for the symmetric eigenproblem. _Linear Algebra and its Applications_, 34, 235–258.

Parlett, B. N. (1998). _The Symmetric Eigenvalue Problem_. Society for Industrial and Applied Mathematics.

Parlett, B. N., & Scott, D. S. (1979). The tracking of spurious eigenvalues in the Lanczos algorithm. _Mathematics of Computation_, 33(145), 217–238.

Polyak, B. T. (2007). Newton's method and its use in optimization. _European Journal of Operational Research_, 181(3), 1086–1096.

Powell, M. J. D. (1978). A fast algorithm for nonlinearly constrained optimization calculations. In G. A. Watson (Ed.), _Numerical Analysis_ (pp. 144–157). Springer.

Saad, Y. (1980). On the rates of convergence of the Lanczos and the block-Lanczos methods. _SIAM Journal on Numerical Analysis_, 17(5), 687–706.

Shanno, D. F. (1970). Conditioning of quasi-Newton methods for function minimization. _Mathematics of Computation_, 24(111), 647–656.

Simon, H. D. (1984). The Lanczos algorithm with partial reorthogonalization. _Mathematics of Computation_, 42(165), 115–142.

Wilkinson, J. H. (1961). Error analysis of direct methods of matrix inversion. _Journal of the ACM_, 8(3), 281–330.

Wilkinson, J. H. (1965). _The Algebraic Eigenvalue Problem_. Oxford University Press.
