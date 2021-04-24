### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# ╔═╡ 17fc4346-a462-11eb-0484-7b4029b01f4d
using Images

# ╔═╡ 3ca7a180-e6ad-4734-b19f-4480de52d13c
md"
# Dataset: Water Clusters as Graphs
"

# ╔═╡ 30339b9c-95a9-4a0e-b7d6-abd65f958e3d
md"
The database we will be working with primarily comes from my own research. In the past, we have generated a very large database of molecular structures called **water clusters**. Water clusters are gas-phase blobs of water molecules which are held together by **hydrogen bonds**. They can consist of essentially any number of water molecules. The actual database consists of the cartesian coordinates of each atom in 3D space and can be accessed at [sites.uw.edu/wdbase](https://sites.uw.edu/wdbase/database-of-water-clusters/). Because of the complexity of analyzing this data, I will be working with a simplified version of this database which treats each water cluster as a graph.
"

# ╔═╡ 5bc0d0b1-0078-4f84-baa0-d0267f8bd7b6
load("./assets/Figure_1_a_W3_25_image_combined.png")

# ╔═╡ 2761102e-78c2-4a4a-a7e2-c972c1bd1f4a
md"
The above image shows the **lowest energy** structures of $\mathrm{(H_2O)_n}$ where $\mathrm{n}$ is the number of water molecules present in the cluster. As you can see, the hydrogen bonds connecting each water molecule can be thought of as defining a graph. In these graphs, each node represents a water molecule and each edge represents a hydorgen bond. We can treat these as either simple graphs (no edge direction) or directed graphs in which each edge is directed along the hydrogen bonds, pointing from a hydrogen atom to an oxygen atom.
"

# ╔═╡ 420184cc-1258-4f98-880e-40dc47914011
load("./assets/cluster_to_graph_illustration_annotated.png")

# ╔═╡ 6b3f36e0-efa5-4010-8baa-758b4a72e75d
md"
The above image illustrates this transformation. On the left is a visualization of each atom in a $\mathrm{(H_2O)_{10}}$ cluster. In the middle is a depiction of a directed graph. The edge directions go from the hydrogen atom (white atoms in the left picture) to an oxygen atom (red atoms). Finally, the right is a simple graph representing the structure. This undirected version will sometimes be referred to as a \"family\" because it defines the connectivity of molecules, but there are many different ways the molecules can be conected by hydrogen bonds, and each of these clusters will have a different energy.
"

# ╔═╡ 78efba80-47b8-4c9c-9435-a9034e7d4053
md"
With that context out of the way, let's describe the actual data I've extracted from these graphs. First, I have enumerated the number of primitive cycles in the simple graphs from size three up to size ten. As an example, the above simple graph contains two 5-cycles (the pentagons) and five 4-cycles (the squares which connect the pentamers). However, we do not count any 8-cycles which would be formed by two neighboring 4-cycles. In this sense, the 4- and 5-cycles are the primitives of this graph.

Additionally, we have information from the directed graphs which describes the connectivity of each molecule in a structure. Namely, molecule 1 in the directed graph above is conventionally described as an AAD molecule. This means it accepts two hydrogen bonds and donates one. Molecule 2 would be called ADD since it accepts one hydrogen bond and donates two.

Finally, we have the total energy of each molecule. The energy is a measure of the stability of each water cluster. Roughly speaking, the lower the energy, the more physically relevant the structure is, as it is much more likely to be present in experiments or the atmosphere, etc.
"

# ╔═╡ 21a82ca1-6b0b-489a-82e8-9a19804829e7
md"
## Questions of Interest:
"

# ╔═╡ 8e6b5906-b80a-45d7-aa8b-b73b670f65da
md"
Some questions I would like to be able to answer are the following:
1. How does the distribution of rings change with cluster size? We should expect 6-membered rings to outpace 3-, 4-, and 5-membered rings as clusters get larger because 6-membered rings are dominant in solid ice.

2. Can we confirm the conventional wisdom that more hydrogen bonds will lead to more stable clusters?

3. Are there any interesting trends when relating different hydrogen bonding configurations (e.g. AAD vs. ADD) to the total energy of a cluster?
"

# ╔═╡ 7c8a60f4-896a-4fdc-9971-1fa48d9df27a
md"
### How are rings distributed with cluster size?
"

# ╔═╡ 9e29427e-4967-4d27-a885-0ac6d0319da8
md"
Let's begin by looking at the distribution of these ring cycles with the number of molecules in each cluster.
"

# ╔═╡ 948c2330-ea9b-4350-a276-ebe940340412
load("assets/rings_vs_cluster_size_lowest_family.png")

# ╔═╡ 47f26447-a892-49cc-926f-ca1db8ee64b9
md"
I have plotted the number of each size cycle against the size of the water cluster. The points represent an average across all graphs of that size, and the shaded area is the standard deviation of each distribution.

We can see that as the clusters get larger, the larger rings becomes much more common. This is to be expected just from the combinatorial nature of the number of paths of a given size.

Typically, chemists think of the molecular structures containing rings in terms of 3-, 4-, 5-, and 6-membered rings. For this reason, I have plotted these distributions in the inset. We can see that at smaller cluster sizes, the smaller rings are more prevalent. However, at 21 molecules, 6-cycles become more common and seem to stay more prevalent than 3-, 4-, and 5-cycles. This is the behavior we would hope to see because bulk ice consists of a lattice of 6-cycles, so we should approach that limit, which seems to be the case.
"

# ╔═╡ d6bcbd4e-d4e6-4e7f-a7eb-41fcb7cac2d7
load("assets/ring_size_correlation.png")

# ╔═╡ 7d873195-eeed-4582-a344-9e7c90eae459
md"
Let's see if these cycle sizes have any strong correlations amongst each other. Above is the pair correlations of the counts of each cycle size from 3 to 10. Roughly speaking, it seems that clusters which are more similar in size are more likely to appear at similar rates. Perhaps this is unsurprising.
"

# ╔═╡ 54f8cd88-f42e-42f2-a62f-8149926a21de
md"
### How are different water bonding arrangements distributed with respect to cluster size?
"

# ╔═╡ 1fe20c54-4c54-4f66-893b-89ff3aef08b4
load("assets/population_of_ring_types_with_cluster_size.png")

# ╔═╡ a5a29ab5-fb59-44d7-9db5-356010f00b23
md"
Let's examine the molecule-level data we have for each graph. As a reminder, AAADD means that a molecule would have three acceptors and two donor molecules. In terms of graphs, this would mean a particular node in the graph has in-degree 3 and out-degree two.

Since there are many more nodes in larger clusters, the counts of each occurence are heavily skewed towards towards the larger water clusters. To a certain extent, this is good because these clusters are better models of liquid water and/or ice. On the other hand, we still care about the smaller clusters, so I have plotted the data on a log-plot on the left and on linear axes on the right. The color bar indicates the number of water molecules present in the cluster. It is commonly said that water is tetrahedrally coordinated. This would correspond to the AADD configuration, which is indeed the most common type of node, and it seems to be getting more common as clusters get larger. Nonetheless, AAD and ADD configurations are still fairly common

Somewhat surprisingly, there is a non-negligible number of AAADD molecules. This is surprising because people usually think of water as only being capable of accepting two hydrogen bonds and donating two hydrogen bonds. Further work needs to be done to understand the importance of these structures.
"

# ╔═╡ b35efb02-713a-4795-a249-84cb4d506abf
md"
### How are hydrogen bonds related to cluster stability?
"

# ╔═╡ 153c2888-1cdf-41d3-bafb-8a1bbeaab03a
load("assets/strength_of_hbond_with_cluster_size.png")

# ╔═╡ 6fe153a6-5262-4baa-b0f1-19cf140c3f0a
md"
When we plot the energy of each cluster against the number of hydrogen bonds (sum of degrees of each node), we get the above plot. First of all, the banding occures because there are discrete jumps in energy when the cluster size increases by one. Since hydrogen bonds are the stabilizing interaction between water molecules, conventional wisdom would say that the more hydrogen bonds are present, the more stable the cluster will be. The above plot seems to show that at lower cluster sizes this is the case, but as the cluster gets bigger, there is a wider range of number of hdyrogen bonds which correspond to low energies. This can be seen by the growing width of each band.

It may be illuminating to normalize these numbers in some way. Let's try that.
"

# ╔═╡ 34627d00-e9fc-4357-b626-4c5c0b1ed6b4
load("assets/strength_of_hbond_with_cluster_size_normalized.png")

# ╔═╡ 7c5a19cd-9952-4ecc-adef-676031112e99
md"
In the above two visualizations, I have normalized the x-axis by the number of molecules in the cluster (this is a common normalization in the literature). The plot on the left normalizes the y-axis by the energy, while the plot on the right normalizes by the number of hydrogen bonds.

These plots are both very illuminating in slightly different ways. The figure on the left shows that the binding energy per molecule increases until there are about three hydrogen bonds per molecule. For reference, liquid water has an average of about 3.6 hydrogen bonds per molecule, and ice has exactly 4 per molecule. Once clusters are large enough to have this many hydrogen bonds per molecule, the binding energy seems to level off. This can be thought of as the clusters reaching some bulk limit.

One might be led to believe that the hydrogen bonds are getting stronger with cluster size by looking at the plot on the left, but in fact, the plot on the right indicates that if anything, there is a slight decrease in the average hdyrogen bond strength with cluster size. That is, just as we have seen above, the number of hydrogen bonds in a cluster is not necessarily a good preictor of the stability of a cluster because accomodating more hydrogen bonds might require adopting sub-optimal configurations and hence result in weaker hydrogen bonds.
"

# ╔═╡ b920bfbb-82ad-4f57-8e1c-f0001e7d77a0
md"
Having looked at the correlation between hydrogen bonds and energy, let's see if there's any similar correlation between the type of molecules in a cluster and that cluster's energy.
"

# ╔═╡ 0eda0aaf-63b7-4540-af84-667ca53d57b0
md"
### How are bonding configurations of water related to cluster stability?
"

# ╔═╡ c9b0f626-65c3-4912-9b66-e76f7decc0f6
load("assets/rings_vs_energy_20_to_30.png")

# ╔═╡ aa2ee694-5f1e-47ac-93e7-f3f99023103f
md"
For visual clarity, we do this analysis for clusters of 20 to 30 molecules. Once again, we get distinct bands corresponding to cluster size. The most obvious trend is that having zero AD molecules is definitely the most favorable situation for finding a low-energy cluster.

Note that the dashed line is the lowest energy structure at each cluster size.

For ADD and AAD molecules, it seems to be that being somewhere in the middle of the distribution is certainly the most favorable. This makes sense given the trade-off between number of hydrogen bonds and overall hydrogen bond strength which we just saw.

Perhaps the most interesting trend is that the number of AADD molecules doesn't seem to matter nearly as much as having an even number of AADD molecules. This is much easier to see if we plot this by itself.
"

# ╔═╡ ef7eb7e2-492e-4203-9106-90342735e47c
load("assets/rings_vs_energy_20_to_30_AADD_only.png")

# ╔═╡ 5e708228-6075-47f2-9153-991e8dcfb3f3
md"
While it's not the most pronounced trend ever, it's now clearer that having 4, 6, 8, or 10 AADD molecules is likely to be more stable than having 2, 3, 5, 7, 9, or 11. Two is likely just too few hydrogen bonds overall. The odd-numbered structures must bias the structures in some way as to make the connectivity less favorable. This is an intriguing but subtle trend which will have to be investigated further in the future.
"

# ╔═╡ e065e882-a5cb-4c5f-bead-434f839d909a
md"
## Conclusions:
"

# ╔═╡ 48d3550d-c965-47c3-995e-d8c6fc69606e
md"
To answer the questions we set out with:
1. We were able to show that 6-cycle rings do indeed outpace 3-, 4-, and 5-cycles as clusters get larger. Even larger rings seem to get more important. Further investigation of these rings is required though because these are not rings which chemists usually think of as being important.

2. We did see the trend that more hydrogen bonds generally correlates with a more stable structure, but we found that as clusters get larger, the number of hydrogen bonds can vary quite a bit while still giving a similar total energy. In the future, identifying the soecific arrangement which are most stable would be quite interesting.

3. Finally, we were able to tease out interesting trends associated with the number of different hydrogen bonding configurations and the stability of the cluster. Namely, AADD configurations seem to be most stabilizing when there are an even number of configurations. This was unexpected and merits further investigation. It was also clear that having AD molecules is generally quite bad for having the lowest energy at a given cluster size.
"

# ╔═╡ 7df79f4b-14e0-483a-afa9-30e66cc4d537
md"
Some comments on data quality:

Fortunately, because this data comes from my own research, I was able to guarantee that there were not any missing values in the dataset. On the other hand, the definition of rings within a cluster is somewhat arbitrary, so the abundance of large ring sizes may be an artifact of the way the analysis is carried out, and hence this may change over time.
"

# ╔═╡ Cell order:
# ╟─17fc4346-a462-11eb-0484-7b4029b01f4d
# ╟─3ca7a180-e6ad-4734-b19f-4480de52d13c
# ╟─30339b9c-95a9-4a0e-b7d6-abd65f958e3d
# ╟─5bc0d0b1-0078-4f84-baa0-d0267f8bd7b6
# ╟─2761102e-78c2-4a4a-a7e2-c972c1bd1f4a
# ╟─420184cc-1258-4f98-880e-40dc47914011
# ╟─6b3f36e0-efa5-4010-8baa-758b4a72e75d
# ╟─78efba80-47b8-4c9c-9435-a9034e7d4053
# ╟─21a82ca1-6b0b-489a-82e8-9a19804829e7
# ╟─8e6b5906-b80a-45d7-aa8b-b73b670f65da
# ╟─7c8a60f4-896a-4fdc-9971-1fa48d9df27a
# ╟─9e29427e-4967-4d27-a885-0ac6d0319da8
# ╟─948c2330-ea9b-4350-a276-ebe940340412
# ╟─47f26447-a892-49cc-926f-ca1db8ee64b9
# ╟─d6bcbd4e-d4e6-4e7f-a7eb-41fcb7cac2d7
# ╟─7d873195-eeed-4582-a344-9e7c90eae459
# ╟─54f8cd88-f42e-42f2-a62f-8149926a21de
# ╟─1fe20c54-4c54-4f66-893b-89ff3aef08b4
# ╟─a5a29ab5-fb59-44d7-9db5-356010f00b23
# ╟─b35efb02-713a-4795-a249-84cb4d506abf
# ╟─153c2888-1cdf-41d3-bafb-8a1bbeaab03a
# ╟─6fe153a6-5262-4baa-b0f1-19cf140c3f0a
# ╟─34627d00-e9fc-4357-b626-4c5c0b1ed6b4
# ╟─7c5a19cd-9952-4ecc-adef-676031112e99
# ╟─b920bfbb-82ad-4f57-8e1c-f0001e7d77a0
# ╟─0eda0aaf-63b7-4540-af84-667ca53d57b0
# ╟─c9b0f626-65c3-4912-9b66-e76f7decc0f6
# ╟─aa2ee694-5f1e-47ac-93e7-f3f99023103f
# ╟─ef7eb7e2-492e-4203-9106-90342735e47c
# ╟─5e708228-6075-47f2-9153-991e8dcfb3f3
# ╟─e065e882-a5cb-4c5f-bead-434f839d909a
# ╟─48d3550d-c965-47c3-995e-d8c6fc69606e
# ╟─7df79f4b-14e0-483a-afa9-30e66cc4d537
