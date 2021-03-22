### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ d20942a4-51f2-11eb-0977-33f2082de974
md"
# TroPrecLS: Tropical Precipitation over Land and Sea

This project investigates the nature of precipitation dynamics in the tropics over land and sea, and if there are any differences in these dynamics to be found in different geographical areas.

I have created several different notebooks to help better illustrate my thought process during the course of this project.  These notebooks will be mostly used for data analysis and plotting.  For computationally heavy, memory intensive, and lengthy work, scripts will be provided and referenced as need be, so that the user can submit these jobs to a backend/cluster to obtain the data themselves.

The general outline of this project is:
1. Defining domains within the tropical region to be explored

2. Download and exploratory analysis of satellte rainfall data in the tropical domain and subdomains
* Analysis of mean, diurnal cycle in domains

3. Download and exploratory analysis of reanalysis ERA5 data in the tropical domain and subdomains
* Analysis of mean, diurnal cycle in domains
* Analysis of vertical profiles for pressure-level data

4. Investigating the relationship between column saturation fraction $r$ and precipitation rate $P$ in reanalysis data
* What is the difference between land and sea
* Is there a difference in this landmass over different land areas

5. Extracting Large-Scale Forcing for SAM
* Quick exploratory analysis of moisture flux convergence
* Finding domain-averaged moisture flux convergence for large-scale forcing input

6. Implementing the Weak-Temperature Gradient Approximation in SAM
* Finding an equilibrium RCE sounding profile
* Transitioning the model from RCE/pseudo-RCE to WTG state
* Investigating the senstivity of WTG momentum damping parameter $a_m$

7. Explore the equilibrium RCE states in SAM
* Plot the 2D rain/PWV fields ... is there convective aggregation?


8. The Impact of the Diurnal Cycle on P-r
* Calculate energy imbalance when slab depth is infinity
* Investigate slab-depth and the diurnal amplitude of SST
"

# ╔═╡ Cell order:
# ╟─d20942a4-51f2-11eb-0977-33f2082de974
