<img src="https://user-images.githubusercontent.com/49268378/87043206-ff710780-c1f4-11ea-81c5-36523124ddf0.jpg" width="550" height="150"> 
TALKIEN - crossTALK bIpartite Network -

You can plot customized networks, obtain network parameters and enrichment analysis.
INPUT

User is required to input two list of genes/proteins. Each of them belongs to one of the two entities the user in interested in. To decipher the molecular crosstalk between them, TALKIEN only uses secreted and membrane receptors proteins. The hypothesis behind is that cellular communication could be partially explained by physical interactions between receptor proteins and their specific ligands. Specifically in paracrine signalling, receptors are activated by ligands proteins secreted by another cell. However, TALKIEN also allow users to visualize pathways downstream receptors in the network (see below): 

<img src="https://user-images.githubusercontent.com/49268378/87044332-ac984f80-c1f6-11ea-9c6b-1c2189a5a01f.jpeg">

These lists may derive from a variety of experiments and could be measured with diverse transcriptomics (such as RNA-seq, microarray or scRNA-seq) and proteomics approaches (mass spectrometry or protein chips among others). For example, co-culture experiments of cancer cells with/without stromal cells, or cells growing under the stimulation of supernatant media collected from immune cells cultures. Lists from the example came from the following transcriptomic experiment: 

<img src="https://user-images.githubusercontent.com/49268378/87044646-216b8980-c1f7-11ea-840a-53570c701504.png" witdh="350" height="350">

PLOT NETWORKS

The basic elements of the networks are nodes (Proteins) and edges (interactions between them) that connect nodes within the network. Connections are defined according to a score based on certain evidence of interaction.
User options.

1. Upload lists Input List #1 and Input List #2. Upload files with header and one column lenght. Otherwise column #1 will be used as reference. Both files must have the same annotation: entrez ID, ensembl ID, gene symbol or uniprot ID are allowed.

<img src="https://user-images.githubusercontent.com/49268378/87044850-68f21580-c1f7-11ea-901f-55e3175934d5.png" width="350" height="150">

There is a preloaded example with two lists. By selecting “Load example data” option, an example network will be displayed. 

<img src="https://user-images.githubusercontent.com/49268378/87045008-9fc82b80-c1f7-11ea-9be9-db4b6ec6a971.png">

2. Choose **Network Database**. By default TALKIEN will combine the two different databases to give wider results.
  * STRING
  * HIPPIE
  * both

<img src="https://user-images.githubusercontent.com/49268378/124778201-169a2080-df41-11eb-8e8c-77cf5d2214b0.png" width="350" height="150">

3. Choose **information source**. Nodes have been classified into plasma membrane **Receptors** or extracellular **Secreted** proteins depending on two different sources. Some divergency may be observed between each other.
  * Human Protein Atlas
  * Uniprot

<img src="https://user-images.githubusercontent.com/49268378/124778890-9e802a80-df41-11eb-99b5-7b938330e2d5.png" width="350" height="150">

4. Choose annotation type.

- Gene Symbol
- Entrez
- Ensembl
- Uniprot

<img src="https://user-images.githubusercontent.com/49268378/87045210-e74eb780-c1f7-11ea-8939-699b63014aee.png" width="175" height="140">

5. Choose network type. Users are allowed to compute interactions between lists (the so-called crosstalk) or interactions between and within lists.

- whole –> get interactions from list 1 to list1, from list2 to list2, from list1 to list2, and viceversa
- crosstalk –> get interactions only from list1 to list2 and viceversa

<img src="https://user-images.githubusercontent.com/49268378/87045423-27ae3580-c1f8-11ea-8030-d618b0b81860.png" width="150" height="75">

 6. Select **score threshold**. That is the minimum interaction value to link two proteins. Only protein-protein interaction with scores higher than chosen will be shown in the network. Score threshold is based either on STRING's combined score or HIPPIE weight score.

<img src="https://user-images.githubusercontent.com/49268378/87045626-6348ff80-c1f8-11ea-9142-63735224f289.png" width="350" heigth="150">

7. Show **downstream interactions**.
  If selected, all interactions coming from the same list of its receptor parent will be shown. This tree-like structure allow users to have an idea of the possible activation cascades switched on by specific receptors

<img src="https://user-images.githubusercontent.com/49268378/124779372-033b8500-df42-11eb-91f0-62e39c1b8a91.png">

  8. Choose **graphical layout**.
Users are allowed to choose between 5 different layouts.
  * Circular layout.
  * Force-directed layout (Fruchterman Reingold algorithm).
  * Sphere layout.
  * Large Graph layout. Preferred option for large networks. -only for network type **whole**-.
  * Bipartite layout. -only for network type **crosstalk**-.

<img src="https://user-images.githubusercontent.com/49268378/87045880-b15e0300-c1f8-11ea-9ea2-10fd71ac9d48.png" width="150" height="125">


All results could be downloaded by clicking donwload buttons on the bottom of the tabs. Network plots could be downloaded first as html objects. Additionally, there is a black button to export the image to a png format once the downloaded html has been opened. 

<img src="https://user-images.githubusercontent.com/49268378/87046925-18c88280-c1fa-11ea-9bbd-f04d84545737.png" width="400" height="350">

MAIN TAB: PLOTS

Network will be displayed and able to be downloaded as an interactive html object. The size of the node is proportional to its degree and the width of the edge is proportional to its interaction score. In **bipartite** mode. There is the option of changing network orientation (UP-DOWN or LEFT-RIGHT) and changing the groups (bipartite network by lists or by location)

<img src="https://user-images.githubusercontent.com/49268378/87046571-a9529300-c1f9-11ea-96f5-be59aed316a9.png" width="250" height="150">

Users can select one node by clicking on it. After selecting one of them, a table with the node's belonging pathways will be displayed. It shows pathway ID, pathway name and a link to see more in-depth pathway properties on reactome website. Moreover, if **downstream interactions** mode is selected, three additional columns will show pathway's size, the number of present nodes in network belonging to each of the displayed pathways and a F-test p-value for the enrichment term.

<img src="https://user-images.githubusercontent.com/49268378/124780601-05521380-df43-11eb-9b4a-b3a98f470345.png">

NETWORK PARAMETERS TAB

For the plotted network, the network descriptive parameters are displayed: nodes, edges, diameter, shortest path, density, average neighbors, clustering coefficient and centralization.

  * nodes: total number of genes/proteins
  * edges: number of links between nodes (total number of interactions)
  * diameter: shortest distance between the two most distal nodes
  * shortest_path: average minimum number of edges between a pair of nodes
  * density: proportion of edges divided by the maximum number of edges the network can have (if all nodes were connected between themselves)
  * average neighbors: average number of links for each node
  * clustering coefficient: number of closed triplets (three nodes completely linked) over the total number of triplets (completely or uncompletely linked)
  * centralization: Average degree centrality


For all nodes, five different centrality methods are displayed and can be downloaded. 

  * Degree: Number of links incident upon a given node
  * Closeness: Average length of the shortest path between a given node and all the other nodes
  * Betweenness: Number of times a given node belongs to a shortest path between two other nodes
  * Eigenvector: Score based on the degree of its direct linked nodes.
  * Page-Rank: Variant of Eigenvector taking into account the direction of links

For each pair of connected nodes, interaction score, location, list and other annotation IDs are displayed and can be downloaded.

ENRICHMENT ANALYSIS TAB

Users are able to perform a functional analysis based on Reactome Pathway Database by selecting this tab

Enrichment analysis results are displayed in two different ways.

  * Dotplot with up to 15 top enriched pathways.
  * Table summarizing enriched pathways, proteins involved, p-values, enrichment scores, ect.
