# TALKIEN - crossTALK InteractivE Network

<p align="justify">
A tool for plot, analize and perform functional analysis molecular crosstalk. The hypothesis behind is that cellular communication could be partially explained by physical interactions between receptor proteins and their specific ligands. In paracrine signalling, receptors are activated by lingand proteins secreted by another cell. TALKIEN also allow users to visualize activated pathways downstream receptors.

[Link to TALKIEN](https://shiny.odap-ico.org/talkien/)

### INPUT

A single tab delimited file with two columns: First column must contain genes or proteins. Second column must contain the tissue, cell type or condition. Columns should be named, otherwise first row will be considered the header. File should look like this:

These lists may derive from a variety of experiments and could be measured with diverse transcriptomics, or proteomics approaches (RNAseq, scRNAseq, Spatial Transcriptomics, microarray, mass spectrometry, among others). For example, co-culture experiments of cancer cells with/without stromal cells, or cells growing under the stimulation of supernatant media collected from immune cell cultures. List from the example came from the following transcriptomic experiment:

<img src="/home/46962313Q/talkien/suppl/fig2_exampleData.jpg", width = 400/>

### NETWORKS

The basic elements of the networks are nodes (proteins) and edges (interactions between them) that connect nodes within the network. Connections are defined according to:

  1. Interaction between a Ligand and a Receptor present in at least one of the Ligand-Receptor DataBases (ddbb)
  2. Interaction between a Receptor and a downstream protein present in STRING ddbb with an interaction score based on certain evidence of interaction

**User options.**

  1. Upload file **Input file**. Files must be uploaded with header and two columns wide with tab delimitaions. First column will contain genes/proteins, second column tissues/cell types/conditions.

There is a preloaded example data with two lists. By selecting "Load example data" option, an example network will be displayed.

  2. When TALKIEN identifies more than two different tissues/lists, all pairwise possibilities are computed by default. However, users can select any combination of pairwise interactions between input lists (up to 8 different lists). Unselection of all possibilities gives the same result than selecting all of them.


  3. Choose **Ligand-Receptor DDBB**. By default, TALKIEN will combine all ligand-receptor DDBBs to integrate all possible information. This will give, in general, wider results.

  4. Choose **annotation type**. Users must provide files with entries in one of the specified annotations:

  5. Choose **network type**. Users can choose between selecting only paracrine interactions (between tissues: crosstalk mode) or both paracrine and autocrine interactions (between and within tissues: full mode)

  6. Compute **downstream interactions**. If selected, all interactions coming from the same list of their parent receptors will be added to the ligand-receptor network. This tree-like structure allow users to have an idea of the possible activation cascades switched on by specific receptors

  6.1 **Curated interactions**. STRING protein-protein interactions (PPI) annotated as experimentally determined or extracted from curated DDBBs. Also an additional filtering based on STRING combined score > 0.3 is done. This downstream network is based on curated interactions.
  6.2 **High confident interactions**. STRING protein-protein interactions (PPI) with STRING combined score > 0.95. According to STRING:"The combined score is computed by combining the probabilities from the different evidence channels and corrected for the probability of randomly observing an interaction.". Thus, by subseting interactions with a combined score > 0.95, we expect high-confident interactions.

  7. Choose **graphical layout**. Five different layouts are available. Large graph layout is the best option for very large connected networks.  Bipartite layout is only available for **crosstalk mode**.

 ***
 

### RESULTS

Output results are organized in five tabs.

##### PLOT TAB
An Interactive graphical network representation will be rendered after choosing the desired options. Users can customize some network elements or filter interactions. Click on **customize plot** to change size of the nodes depending on centrality measures or showing only specific components of the network whenever possible.

Searching nodes is possible either by selecting from an id list or by mouse hover over each node

Node manipulation

 * Mouse hover over a node: highlight of its name
 * one-click over a node (when downstream analysis is switched off): highlight its direct connections
 * one-click over a node (when downstream analysis is switched on): highlight its direct connections (for ligands), or highlight all direct and indirect connections between a receptor and downstream nodes (for receptors and downstream nodes)
 * Mouse hold button over a node: allows to manually move the node to another position.

After node selection, all of the Reactome pathway IDs in which the node is involved appear. There are links to Reactome website for more details. If downstream analysis is switched on, a p-value is computed taking into account genes in the input lists. In this scenario, only Reactome pathway IDs with an adjusted p-value smaller than 0.05 will be printed. In addition, the length of the pathway and the highlighted nodes present in the pathway will also be printed.

***
  
##### NETWORK INTERACTIONS TAB
For each pair of connected nodes, tissue/cell type of origin, and protein type (ligand, receptor, or downstream) are displayed.

***
 
##### NETWORK PARAMETERS TAB
For the plotted network, descriptive and topological parameters of the network are displayed:

 * Nodes: number of nodes in the network
 * Edges: number of connections between all nodes in the network
 
 next parameters will be computed for giant component if the network is not fully conected:
 
 * Diameter: maximum distance between the further nodes in the network
 * Shortest path: average minimum distance between two different nodes
 * Density: number of connections relative to the total possible connections in the network
 * Average neighbors: average number of connections for each node
 * Clustering coefficient: measures the degree to which nodes in a network tends to cluster together
 * Centralization: average normalized degree of the network
 * Components: number of isolated networks.


For each node, centrality measures, tissue/cell type of origin, protein type (ligand, receptor, or downstream) and all different annotation IDs are displayed and can be downloaded. Centrality measures computed in TALKIEN are:

 * Degree: number of edges incident upon a node
 * Closeness: Average length of the shortest path between a node and all other nodes in the graph
 * Betweenness: number of times that a node is present along the shortest path between two nodes
 * Eigenvector: Measures the influence of nodes based on the importance of its direct neighbors
 * Page-Rank: Slightly different centrality measure than eigenvector for directed networks
 * Eccentricity: maximum distance from a node to any other node in the network

Other network measures are computed for each node:

 * Clustering coefficient: measures the density of the connections between the direct neighbors of a node
 * Component: the number of the component in which a node is present

***

##### ENRICHMENT ANALYSIS TAB
Users are able to perform functional analysis based on Reactome Pathway Database by selecting this tab. There are some extra options to customize the analysis.

 * When downstream analysis is not selected enrichment is done against ligands and receptors.
 * When downstream analysis is selected, enrichment is done against all annotated entries
 * In addition, is also possible to perform enrichments for specific network clusters, tissues/lists

Enrichment analysis results are displayed in two different ways.

 * Dot plot with top enriched pathways.
 * Table summarizing enriched pathways: pathway ID, pathway name, gene ratio (number of selected nodes relative to the number of all the nodes in the network), background ratio (number of pathway elements relative to the background), p-value, adjusted p-value, gene names for the nodes in the network present in the pathway, and links to specific Reactome Pathways.

***

##### DRUGS TAB
Network receptors with available knowledge of drug interactions (drug-gene interactions) from DGIdb are displayed. Table contains different annotations for each receptor, name of the drug and alternative ids, interaction scores between gene and drug, source of the interaction and PMIDs where the interaction has been reported.

***
All results can be downloaded.

Code freely available at:
https://github.com/odap-ubs/talkien/
</p>
