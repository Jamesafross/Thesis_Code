using Plots

τ = 1
N = 1000
R = LinRange(0,5,N)
V = LinRange(-10,10,N)

sync = zeros(N,N)
angle1= zeros(N,N)
for i = 1:N
    for j = 1:N
        W = τ.*pi.*R[i] .+ im.*V[j]
        Z = (1 .-conj.(W)) ./(1 .+conj.(W))
        sync[i,j] = abs.(Z)
        angle1[i,j] = angle.(Z)
    end
end

p1 = heatmap(sync')
p2 = heatmap(mod.(angle1',2pi))
plot(p1,p2)
