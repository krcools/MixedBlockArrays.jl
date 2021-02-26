using MixedBlockArrays
using SparseArrays
using LinearAlgebra

ax1 = blockedrange([1,2])
ax2 = blockedrange([1,1,2])
T = Float64
A = MixedBlockArray{T}((ax1,ax2))

A[Block(1,3)] = Fill(2, (1,2))
A[Block(2),Block(2)] = SparseMatrixCSC{Float64}(I, 2, 1)
A[Block(2,2)] = zeros(2,2) 

A[Block(2,2)][1,1] = pi
A[3, 2] = -pi

@assert A == [0 0 2 2; 0 pi 0 0; 0 -pi 0 0]