# Toolkit for simulating 2D dynamical system on surface
module Dynamical2D

export D, ∇, divF, curlF, diffOpsM, kutta, SurfaceFun, SetShape

using LinearAlgebra, ForwardDiff

module SurfaceFun
    # Quadric surface of the form z = -(1/r)*[x,y]'*A*[x,y]
    quadric(x,y;r=4.0,A=[1 0; 0 1]) = [x,y,-(1/r)*[x,y]'*A*[x,y]]
    # Shpere of θ: Longitude, ϕ: Latitude (geographic coodinate)
    sphere(θ,ϕ;r=1.0) = r*[cos(θ)*cos(ϕ),sin(θ)*cos(ϕ),sin(ϕ)]
    # Torus of θ: Longitude (lateral), ϕ: Meridian (vertical)
    torus(θ,ϕ;r=0.5,R=1.0) = [(R+r*cos(ϕ))*cos(θ),(R+r*cos(ϕ))*sin(θ),r*sin(ϕ)]
    # Cone
    cone(x,y;r=1,h=1,c=0) = [x,y,h*(r-sqrt(x^2+y^2))/r+c]
end # module Surface

struct SetShape
    type::Symbol    # should matche to above function names
    f::Function     # function of shape
    par::NamedTuple # parameters
    SetShape(type;par...) = new(
        type,
        (u,v) -> getfield(SurfaceFun,type)(u,v;par...),
        NamedTuple(par)
    )
end

# Projection on torus
function ontorus(p;r=0.5,R=1.0)
    Rvec = R*[normalize(p[1:2]); 0.0]
    rvec = r*normalize(p - Rvec)
    Rvec + rvec
end

# helper function to split args of f
spa(f::Function) = x -> f(x...)

# -π/2 rotation
const Q = [0.0 1.0; -1.0 0.0]

# Differential operators on flat space: Function -> Function
∇(f::Function) = (x,y) -> ForwardDiff.gradient(spa(f),[x,y])
D(f::Function) = (x,y) -> ForwardDiff.jacobian(spa(f),[x,y])
divF(f::Function) = (x,y) -> tr(D(f)(x,y))
curlF(f::Function) = (x,y) -> tr(Q*D(f)(x,y)) # = ∂₁f₂ - ∂₂f₁
#curl(f::Function) = (x,y) -> let J=D(f)(x,y); J[2,1]-J[1,2] end

# Generate differential operators on surface M
function diffOpsM(shapefun::Function; drgi=nothing)
    # Jacobian matrix of R² (param.) → R³ (surface)
    J(x,y) = ForwardDiff.jacobian(spa(shapefun),[x,y])
    # Riemannian metric
    G(x,y) = let J = J(x,y); J'J end
    # Manual definition of |1/√g| is recommended for speed
    # e.g., r^2*cos(y) for the sphere of radius r
    if isnothing(drgi)
        drgi = (x,y) -> 1/det(G(x,y)^(1/2))
    end
    # Derivatives w.r.t. the metric: D(x,y)[:,:,1][i,j] = ∂gᵢⱼ/∂x
    Dg(x,y) = reshape.(eachcol(ForwardDiff.jacobian(spa(G),[x,y])),2,2)
    # Christoffel symbols of the second kind
    # returns (Γ¹ᵢⱼ, Γ²ᵢⱼ) = (g¹ᵐ[ij,m], g²ᵐ[ij,m])
    function Γ(i,j)
        (x,y) -> let ∂g = Dg(x,y)
            G(x,y)\(0.5*[∂g[j][i,m]+∂g[i][j,m]-∂g[m][i,j] for m in 1:2])
        end
    end
    # Covariant derivative of f: usual deriv. + correction 
    Dc(f) = (x,y) -> D(f)(x,y) - [f(x,y)'Γ(i,j)(x,y) for i in 1:2, j in 1:2]
    # Vector calculus on M
    divM(f) = (x,y) -> tr(G(x,y)\Dc(f)(x,y)) # trace w.r.t. the metric
    curlM(f) = (x,y) -> drgi(x,y)*curlF(f)(x,y) # area-scaled curl
    #det(rtgi(x,y))*tr(Q*D(fp)(x,y))
    (J=J,G=G,Dc=Dc,divM=divM,curlM=curlM)
end

# Lunge-Kutta 4: returns s+1 point of the trajectory
function kutta(f::Function, x0::Vector{Float64}; s=100::Int, dt=0.01::Float64)
    X = zeros(Float64,length(x0),s+1) # +1 is required for t=0
    X[:,1] = x0
    for i in 1:s
        k1 = f(X[:,i])dt
        k2 = f(X[:,i] + k1/2)dt
        k3 = f(X[:,i] + k2/2)dt
        k4 = f(X[:,i] + k3)dt
        X[:,i+1] = X[:,i] + (k1 + 2k2 + 2k3 + k4)/6
    end
    return X
end

# Relocation and splitting trajectory in periodic domain
function reformtraj(traj)
    traj = rem2pi.(traj,RoundNearest)
    ep = findall(norm.(eachcol(diff(traj,dims=2))) .> 1)
    traj = Array.(eachcol(traj))
    for p in ep
        insert!(traj,p+1,[NaN,NaN])
        ep .= ep.+1
    end
    traj
end

end # module Dynamical2D