module MixedBlockArrays


using Reexport

@reexport using BlockArrays
@reexport using FillArrays

import Base.Cartesian: @nloops, @ntuple
import BlockArrays: block, blockindex

export MixedBlockArray

struct MixedBlockArray{T, N, BS<:NTuple{N,AbstractUnitRange{Int}}} <: AbstractBlockArray{T, N}
    blocks::AbstractArray{AbstractArray{T,N},N}
    axes::BS

    MixedBlockArray(blocks::R, blocked_axes::BS) where {T, N, R<:AbstractArray{AbstractArray{T,N},N}, BS<:NTuple{N,AbstractUnitRange{Int}}} =
        new{T, N, BS}(blocks, blocked_axes)
end

@generated function MixedBlockArray{T}(axes::NTuple{N,AbstractUnitRange{Int}}) where {T, N, R<:AbstractArray{AbstractArray{T,N},N}}
    return quote
        n_blocks = map(blocklength, axes)
        blocks = Array{AbstractArray{T,N},N}(undef, n_blocks)
        @nloops $N block dim->blockaxes(axes[dim],1) begin
            block_indices = @ntuple $N block
            indices = getindex.(axes, block_indices)
            block_size = length.(indices)
            blocks[Int.(block_indices)...] = Zeros{T}(block_size...)
        end
        return MixedBlockArray(blocks, axes)
    end
end

BlockArrays.getindex(A::MixedBlockArray, b1::Block{1}, b2::Block{1}, btail::Vararg{Block{1}}) = getindex(A, Block(b1.n..., b2.n...), btail...)
BlockArrays.getindex(A::MixedBlockArray, I::Block) = A.blocks[I.n...]

function BlockArrays.setindex!(A::MixedBlockArray, v,b1::Block{1}, b2::Block{1},
    btail::Vararg{Block{1}})

    setindex!(A, v, Block(b1.n..., b2.n...), btail...)
end
BlockArrays.setindex!(A::MixedBlockArray, v, I::Block) = setindex!(A.blocks, v, I.n...)

Base.axes(A::MixedBlockArray) = A.axes
Base.size(A::MixedBlockArray) = length.(axes(A))
function Base.getindex(A::MixedBlockArray, I::Vararg{Int})
    Bi = findblockindex.(axes(A), I)
    B, i = block.(Bi), blockindex.(Bi)
    bl = A[B...]
    return bl[i...]
end

function Base.setindex!(A::MixedBlockArray, v, I::Vararg{Int})
    Bi = findblockindex.(axes(A), I)
    B, i = block.(Bi), blockindex.(Bi)
    bl = A[B...]
    setindex!(bl, v, i...)
end


end
