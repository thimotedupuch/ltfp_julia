### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ d447f068-2892-11f1-adc0-1b0394365f8b
using Pkg; Pkg.activate(".."); Pkg.status()

# ╔═╡ e34a67d0-8c8c-4672-8712-a013298846eb
using CairoMakie, Statistics

# ╔═╡ 188b111c-0a94-4e73-ae03-199af8d51fee
CairoMakie.set_theme!(theme_latexfonts())

# ╔═╡ e937fe95-dd85-4b36-b29b-1752f5f05d18
CairoMakie.activate!(type = "svg")

# ╔═╡ f5610559-4903-4437-8c06-ea5f88b6ca4c
colors = Makie.wong_colors()

# ╔═╡ 492b3db3-e4cc-4a6a-bac3-48608354f5b3
blue = colors[1]

# ╔═╡ 10674023-6d55-4f53-8bfe-d9334c63afa5
orange = colors[2]

# ╔═╡ 567b47b9-e2ef-4494-ad57-c548a1f49518
red = colors[6]

# ╔═╡ 54042c65-3dd6-4310-9bba-3ca4d1e654f8
cummax(x::AbstractVector) = accumulate(max, x)

# ╔═╡ b4ffa839-7dd2-477d-8a98-1b829d727ea4
cummax(A; dims=1) = mapslices(x->accumulate(max, x), A; dims=dims)

# ╔═╡ 94fcc55f-d58d-4d7b-8ef2-52dd41e9986c
function plot_single_draw(ax, n, nrep)
    X = randn(n, nrep)
    scatter!(ax, X[:, 1], markersize=2, color=blue, label="Sample")
    lines!(ax, cummax(X[:, 1]), color=red, label="Cumulative max")
    ax.xlabel = L"n"
    ax.ylabel = "Single draw"
	axislegend(ax, position=:rb)
end

# ╔═╡ 4ad66892-9b62-4da0-af4e-6a3538876b99
function plot_sample_paths(ax, n, nrep)
    X = randn(n, nrep)
    xvals = sqrt.(log.(1:n))
    for j in 1:nrep
        lines!(ax, xvals, cummax(X[:, j]), color=blue)
    end
    ax.xlabel = L"(\log\ n)^{1/2}"
    ax.ylabel = "Samples of max"
	ax.yticks = -4:4
    return X  # reused by the third panel
end

# ╔═╡ f7ac71b4-0dcb-4c35-bfb5-8d56fb78dbda
function plot_expected_max(ax, X)
    xvals = sqrt.(log.(1:size(X, 1)))
    cm = cummax(X, dims=1)
    μ = vec(mean(cm, dims=2))
    σ = vec(std(cm, dims=2))	
    lines!(ax, xvals, μ)
    lines!(ax, xvals, μ .+ σ, linestyle=:dot, color=blue)
    lines!(ax, xvals, μ .- σ, linestyle=:dot, color=blue)
    ax.xlabel = L"(\log\ n)^{1/2}"
    ax.ylabel = "Expectation of max"
	ax.yticks = -4:4
end

# ╔═╡ 6cb7c78e-5315-4172-9897-bea94868683f
function make_figure()
    fig = Figure(size=(900, 300))
	ax1 = Axis(fig[1, 1])
    ax2 = Axis(fig[1, 2])
    ax3 = Axis(fig[1, 3])

    plot_single_draw(ax1, 10_000, 1000)
    X = plot_sample_paths(ax2, 10_000, 10)
    plot_expected_max(ax3, X)

	#save("expectation_of_max_highres.png", fig, dpi=300)
	#save("expectation_of_max_highres.svg", fig)
    
	return fig
end

# ╔═╡ b29b98b0-61a6-468f-b5e4-5ac863c00c06
fig = make_figure()

# ╔═╡ Cell order:
# ╠═d447f068-2892-11f1-adc0-1b0394365f8b
# ╠═e34a67d0-8c8c-4672-8712-a013298846eb
# ╠═188b111c-0a94-4e73-ae03-199af8d51fee
# ╠═e937fe95-dd85-4b36-b29b-1752f5f05d18
# ╠═f5610559-4903-4437-8c06-ea5f88b6ca4c
# ╠═492b3db3-e4cc-4a6a-bac3-48608354f5b3
# ╠═10674023-6d55-4f53-8bfe-d9334c63afa5
# ╠═567b47b9-e2ef-4494-ad57-c548a1f49518
# ╠═54042c65-3dd6-4310-9bba-3ca4d1e654f8
# ╠═b4ffa839-7dd2-477d-8a98-1b829d727ea4
# ╠═94fcc55f-d58d-4d7b-8ef2-52dd41e9986c
# ╠═4ad66892-9b62-4da0-af4e-6a3538876b99
# ╠═f7ac71b4-0dcb-4c35-bfb5-8d56fb78dbda
# ╠═6cb7c78e-5315-4172-9897-bea94868683f
# ╠═b29b98b0-61a6-468f-b5e4-5ac863c00c06
