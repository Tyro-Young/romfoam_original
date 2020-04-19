.. _Theory:

Theory
======

Linear ROM
----------

Linear ROM consists of offline and online stages

Offline stage
~~~~~~~~~~~~~~

- We first run CFD samples with a range of design variable (:math:`x`) values, e.g., varying ramp angle for the Ahmed body

- Then, we assemble state variables from all these CFD samples (:math:`w`; e.g., velocity, pressure, turbulence variables) into a so-called snapshot matrix :math:`S \in \mathbb{R}^{N\times n}`, where :math:`N` is the number of states and :math:`n` is the number of CFD samples, and :math:`N \gg n`. 
 :math:`S^{N \times n}=
 \begin{bmatrix}
 w_{1,1} & w_{1,2} & \cdots & w_{1,n}\\
 w_{2,1} & w_{2,2} & \cdots & w_{2,n}\\
 \vdots & \vdots & \ddots & \vdots \\
 w_{N,1} & w_{N,2} & \cdots & w_{N,n}\\
 \end{bmatrix}`

- Based on the snapshot matrix :math:`S\in \mathbb{R}^{N\times n}`, we compute a set of basis vectors :math:`\Phi \in \mathbb{R}^{N\times n}=[\overrightarrow{\Phi}_1, \overrightarrow{\Phi}_2, \cdots, \overrightarrow{\Phi}_n]` using the Singular Value Decomposition (SVD), where :math:`\overrightarrow{\Phi} \in \mathbb{R}^{N\times 1}` is a basis vector.
 :math:`\Phi^{N \times n}=
 \begin{bmatrix}
 \phi_{1,1} & \phi_{1,2} & \cdots & \phi_{1,n}\\
 \phi_{2,1} & \phi_{2,2} & \cdots & \phi_{2,n}\\
 \vdots & \vdots & \ddots & \vdots \\
 \phi_{N,1} & \phi_{N,2} & \cdots & \phi_{N,n}\\
 \end{bmatrix}`

- With this treatment, a state variable vector $\overrightarrow{w}$ in the design space can be approximated as a linear combination of the basis vectors, i.e., :math:`\overrightarrow{w} \approx \Phi \overrightarrow{w}_r=[w_r^1 \overrightarrow{\Phi}_1, w_r^2 \overrightarrow{\Phi}_2, \cdots, w_r^n \overrightarrow{\Phi}_n]`. The goal of the online stage ROM is to quickly find the reduced vector :math:`\overrightarrow{w}_r` for any given design variable :math:`\overrightarrow{x}`.

Online stage
~~~~~~~~~~~~~

In the online stage, we want to quickly find the reduced vector :math:`\overrightarrow{w}_r` for any given design variable :math:`\overrightarrow{x}`. Then we use :math:`\overrightarrow{w}=\Phi \overrightarrow{w}_r` to get the approximated full states. To this end:

- We first linearize the full-scale nonlinear residuals at a reference point (:math:`\overrightarrow{w}_0`, :math:`\overrightarrow{x}_0`) where we have converged the CFD: 
 :math:`\overrightarrow{R}(\overrightarrow{w}_0+\Delta \overrightarrow{w},\overrightarrow{x}_0+\Delta \overrightarrow{x})=\overrightarrow{R}(\overrightarrow{x}_0,\overrightarrow{w}_0)+\frac{\partial R}{\partial w}_{(\overrightarrow{w}_0,\overrightarrow{x}_0)} \Delta \overrightarrow{w} + \frac{\partial R}{\partial x}_{(\overrightarrow{w}_0,\overrightarrow{x}_0)} \Delta \overrightarrow{x}`

- Since (:math:`\overrightarrow{w}_0`, :math:`\overrightarrow{x}_0`) is a converged solution, :math:`\overrightarrow{R}(\overrightarrow{x}_0,\overrightarrow{w}_0)=0`.

- Now, for a new design variable :math:`\overrightarrow{x}_0+\Delta \overrightarrow{x}`, we need to solve for :math:`\Delta \overrightarrow{w}` such that the residuals at the new point :math:`(\overrightarrow{w}_0+\Delta \overrightarrow{w},\overrightarrow{x}_0+\Delta \overrightarrow{x})` are also zero, i.e., :math:`\overrightarrow{R}(\overrightarrow{w}_0+\Delta \overrightarrow{w},\overrightarrow{x}_0+\Delta \overrightarrow{x})=0`

- In other words, we need to solve
 :math:`\frac{\partial R}{\partial w}_{(\overrightarrow{w}_0,\overrightarrow{x}_0)} \Delta \overrightarrow{w} + \frac{\partial R}{\partial x}_{(\overrightarrow{w}_0,\overrightarrow{x}_0)} \Delta \overrightarrow{x}=0`, or \\ 
 :math:`\frac{\partial R}{\partial w}_{(\overrightarrow{w}_0,\overrightarrow{x}_0)} \Delta \overrightarrow{w} = - \frac{\partial R}{\partial x}_{(\overrightarrow{w}_0,\overrightarrow{x}_0)} \Delta \overrightarrow{x}`

- Let us take a close look at this equation (:math:`n_x` is the number of design variables): 
 :math:`\underbrace{ \frac{\partial R}{\partial w}_{(\overrightarrow{w}_0,\overrightarrow{x}_0)} }_{N \times N} \underbrace{ \vphantom{\dfrac{nan}{nan}} \Delta \overrightarrow{w}}_{N \times 1} = - \underbrace{ \frac{\partial R}{\partial x}_{(\overrightarrow{w}_0,\overrightarrow{x}_0)} }_{N \times n_x} \underbrace{ \vphantom{\dfrac{nan}{nan}} \Delta \overrightarrow{x} }_{n_x \times 1}`

- To reduce this equation, we left-multiple :math:`\Phi^T` to its lhs and rhs and use the relationship :math:`\Delta \overrightarrow{w} = \Phi \Delta \overrightarrow{w}_r`, we obtain:
 :math:`\underbrace{\vphantom{\dfrac{nan}{nan}} \Phi^T}_{n \times N} \underbrace{ \frac{\partial R}{\partial w}_{(\overrightarrow{w}_0,\overrightarrow{x}_0)} }_{N \times N} \underbrace{\vphantom{\dfrac{nan}{nan}} \Phi}_{N \times n} \underbrace{ \vphantom{\dfrac{nan}{nan}} \Delta \overrightarrow{w}_r}_{n \times 1} = - \underbrace{\vphantom{\dfrac{nan}{nan}} \Phi^T}_{n \times N} \underbrace{ \frac{\partial R}{\partial x}_{(\overrightarrow{w}_0,\overrightarrow{x}_0)} }_{N \times n_x} \underbrace{ \vphantom{\dfrac{nan}{nan}} \Delta \overrightarrow{x} }_{n_x \times 1}` 
   
 :math:`\underbrace{ \Phi^T \frac{\partial R}{\partial w}_{(\overrightarrow{w}_0,\overrightarrow{x}_0)} \Phi }_{ A_r \in n \times n} \underbrace{ \vphantom{\dfrac{nan}{nan}} \Delta \overrightarrow{w}_r}_{n \times 1} = - \underbrace{ \Phi^T \frac{\partial R}{\partial x}_{(\overrightarrow{w}_0,\overrightarrow{x}_0)} }_{ B_r \in n \times n_x} \underbrace{ \vphantom{\dfrac{nan}{nan}} \Delta \overrightarrow{x} }_{n_x \times 1}`

 :math:`\underbrace{ A_r }_{n \times n} \underbrace{ \Delta \overrightarrow{w}_r}_{n \times 1} = - \underbrace{ B_r }_{n \times n_x} \underbrace{ \Delta \overrightarrow{x} }_{n_x \times 1}`:  this is an easy :math:`n \times n` linear equation

Summary
~~~~~~~

To summarize, the linear ROM consists of the following steps:

Offline: 

- Generate CFD samples with a range of :math:`x` values (different shapes)
- Perform SVD to compute the :math:`\Phi, A_r`, and :math:`B_r` matrices and save them to the disk

Online: 

- Read in the :math:`\Phi, A_r`, and :math:`B_r` matrices and compute the reduced matrix :math:`A_r` and :math:`B_r` at reference point (:math:`\overrightarrow{w}_0, \overrightarrow{x}_0`)
- For a given new design with :math:`\Delta \overrightarrow{x}`, solve for :math:`\Delta \overrightarrow{w_r}` using :math:`A_r \Delta \overrightarrow{w_r} = -B_r \Delta \overrightarrow{x}`
- Compute the new, full-scale state variables :math:`\overrightarrow{w}` using :math:`\overrightarrow{w}=\overrightarrow{w}_0+\Delta \overrightarrow{w} = \overrightarrow{w}_0 + \Phi \Delta \overrightarrow{w}_r`
- Deform the geometry to the new design point (:math:`\overrightarrow{x}_0+\Delta \overrightarrow{x}`) and compute the new objective functions (e.g., :math:`C_D`) using the new state variables :math:`\overrightarrow{w}`


Nonlinear ROM
-------------

Nonlinear ROM also consists of offline and online stages. 
  
Offline stage
~~~~~~~~~~~~~

The offline stage is similar to the linear ROM.

Online stage
~~~~~~~~~~~~

For the online stage:

- We want to drive the full-scale residual to zero, i.e., :math:`\overrightarrow{R}(\overrightarrow{w},\overrightarrow{x})=0`
- To reduce the above residual, we replace :math:`\overrightarrow{w}` with :math:`\Phi \overrightarrow{w}_r` and left-multiply the residual equation with :math:`\Phi^T`. We obtain: 
 :math:`\underbrace{ \Phi^T}_{n \times N} \underbrace{ \overrightarrow{R}(\Phi\overrightarrow{w}_r,\overrightarrow{x})}_{N \times 1}=\underbrace{\overrightarrow{R}_r(\Phi\overrightarrow{w}_r,\overrightarrow{x})}_{n\times 1}=0`
- Essentially, the reduced nonlinear residual vector :math:`\overrightarrow{R}_r` has size of :math:`n \times 1` and its input variables are :math:`\overrightarrow{w}_r \in \mathbb{R}^{n \times 1}` and :math:`\overrightarrow{x} \in \mathbb{R}^{n_x \times 1}`
 We use a matrix-free Newton-Krylov approach to solve the reduced residual equation :math:`\overrightarrow{R}_r(\Phi\overrightarrow{w}_r,\overrightarrow{x})=0`
- The Newton method says  
 :math:`\frac{\partial \overrightarrow{R}_r}{\partial \overrightarrow{w}_r}_{(\overrightarrow{w}_r^0,\overrightarrow{x}^0)} \Delta \overrightarrow{w}_r^0 = -{\overrightarrow{R}_r}_{(\overrightarrow{w}_r^0,\overrightarrow{x}^0)}`: solved linear euqation with GMRES
 :math:`\overrightarrow{w}_r^1=\overrightarrow{w}_r^0+\Delta \overrightarrow{w}_r^0`: line search to update :math:`\overrightarrow{w}`
 :math:`\frac{\partial \overrightarrow{R}_r}{\partial \overrightarrow{w}_r}_{(\overrightarrow{w}_r^1,\overrightarrow{x}^1)} \Delta \overrightarrow{w}_r^1 = -{\overrightarrow{R}_r}_{(\overrightarrow{w}_r^1,\overrightarrow{x}^1)}` : solved another linear euqation with GMRES 
 .... : repeat until converged

- To enable matrix-free approach, the matrix-vector product :math:`\frac{\partial \overrightarrow{R}_r}{\partial \overrightarrow{w}_r}\overrightarrow{\psi}` used in the above GMRES linear equation solution is computed as   
 :math:`\frac{\partial \overrightarrow{R}_r}{\partial \overrightarrow{w}_r}\overrightarrow{\psi}=\frac{\overrightarrow{R}_r(\overrightarrow{w}_r+\varepsilon \overrightarrow{\psi} )-\overrightarrow{R}_r(\overrightarrow{w}_r)}{\varepsilon}`
 So we don't need to explicitly form and store the :math:`\frac{\partial \overrightarrow{R}_r}{\partial \overrightarrow{w}_r}` matrix, which saves memory and improves speed. Check this paper for more details:  Knoll and Keyes. Jacobian-free Newton-Krylov Methods: a survey of approaches and applications, JCP, 2004.

- The matrix-free Newton-Krylov method does need a preconditioner matrix :math:`\frac{\partial \overrightarrow{R}_r}{\partial \overrightarrow{w}_r}_{\textrm{PC}}` for each GMRES linear equation solution. This PC matrix is computed in the offline stage using the state variables from a reference sample point. It will NOT be updated during the online stage.

- Once we converge the reduced residual :math:`\overrightarrow{R}_r(\Phi\overrightarrow{w}_r,\overrightarrow{x})=0`, we need to also check if the full-scale residual :math:`\overrightarrow{R}(\overrightarrow{w},\overrightarrow{x})=0` also drops (a lot of time they don't!).

- Finally, the full-scale state variables are computed as :math:`\overrightarrow{w}=\Phi \overrightarrow{w}_r`.

Summary
~~~~~~~

To summarize, the nonlinear ROM consists of the following steps:

Offline: 

- Generate CFD samples with a range of :math:`x` values (different shapes)
- Perform SVD to compute the :math:`\Phi` matrix and save it to the disk
- Compute the preconditioner matrix :math:`\frac{\partial \overrightarrow{R}_r}{\partial \overrightarrow{w}_r}_{\textrm{PC}}` and save it to the disk

Online: 

- For a given new design with :math:`\Delta \overrightarrow{x}`, deform the geometry to a new design point (:math:`\overrightarrow{x}_0+\Delta \overrightarrow{x}`)
- Copy the state variable vector from a reference sample  (:math:`\overrightarrow{w}_0`) and use it as the initial condition for the Newton-Krylov solution
- Solve the reduced residual :math:`\overrightarrow{R}_r(\Phi\overrightarrow{w}_r,\overrightarrow{x}_0+\Delta \overrightarrow{x})=0`  using the Newton-Krylov method, and get :math:`\overrightarrow{w}_r`
- Compute the new, full-scale state variables :math:`\overrightarrow{w}= \Phi \overrightarrow{w}_r`
- Compute the new objective functions (e.g., :math:`C_D`) using the new state variables :math:`\overrightarrow{w}`

