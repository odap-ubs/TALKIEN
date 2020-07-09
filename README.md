<img src="https://user-images.githubusercontent.com/49268378/87043206-ff710780-c1f4-11ea-81c5-36523124ddf0.jpg" width="550" height="150"> 
TALKIEN - crossTALK bIpartite Network -

You can plot customized networks, obtain network parameters and enrichment analysis.
INPUT

User is required to input two list of genes/proteins. Each of them belongs to one of the two entities the user in interested in. To decipher the molecular crosstalk between them, TALKIEN only uses secreted and membrane receptors proteins. The hypothesis behind is that cellular communication could be partially explained by physical interactions between receptor proteins and their specific ligands. Specifically in paracrine signalling, receptors are activated by ligands proteins secreted by another cell. However, TALKIEN also allow users to visualize pathways downstream receptors in the network (see below): 

<img src="https://user-images.githubusercontent.com/49268378/87044332-ac984f80-c1f6-11ea-9c6b-1c2189a5a01f.jpeg" width="450" height="150">

These lists may derive from a variety of experiments and could be measured with diverse transcriptomics (such as RNA-seq, microarray or scRNA-seq) and proteomics approaches (mass spectrometry or protein chips among others). For example, co-culture experiments of cancer cells with/without stromal cells, or cells growing under the stimulation of supernatant media collected from immune cells cultures. Lists from the example came from the following transcriptomic experiment: 

<img src="https://user-images.githubusercontent.com/49268378/87044646-216b8980-c1f7-11ea-840a-53570c701504.png" witdh="350" height="350">

PLOT NETWORKS

The basic elements of the networks are nodes (Proteins) and edges (interactions between them) that connect nodes within the network. Connections are defined according to a score based on certain evidence of interaction.
User options.

1. Upload lists Input List #1 and Input List #2. Upload files with header and one column lenght. Otherwise column #1 will be used as reference. Both files must have the same annotation: entrez ID, ensembl ID, gene symbol or uniprot ID are allowed.

<img src="https://user-images.githubusercontent.com/49268378/87044850-68f21580-c1f7-11ea-901f-55e3175934d5.png" width="350" height="150">

There is a preloaded example with two lists. By selecting “Load example data” option, an example network will be displayed. 

<img src="https://user-images.githubusercontent.com/49268378/87045008-9fc82b80-c1f7-11ea-9be9-db4b6ec6a971.png" width="150" height="60">

2. Choose annotation type.

- Gene Symbol
- Entrez
- Ensembl
- Uniprot

<img src="https://user-images.githubusercontent.com/49268378/87045210-e74eb780-c1f7-11ea-8939-699b63014aee.png" width="175" height="140">

Press load data

After clicking load data, results will be displayed.

<img src="https://user-images.githubusercontent.com/49268378/87046490-8f18b500-c1f9-11ea-86fb-ad869e71346f.png" width="150" height="60">

3. Choose network type. Users are allowed to compute interactions between lists (the so-called crosstalk) or interactions between and within lists.

- whole –> get interactions from list 1 to list1, from list2 to list2, from list1 to list2, and viceversa
- crosstalk –> get interactions only from list1 to list2 and viceversa

<img src="https://user-images.githubusercontent.com/49268378/87045423-27ae3580-c1f8-11ea-8030-d618b0b81860.png" width="150" height="75">

4. Select score threshold. That is the minimum interaction value to link two proteins. Only protein-protein interaction with scores higher than chosen will be shown in the network. Score threshold is based on STRING’s combined score. 

<img src="https://user-images.githubusercontent.com/49268378/87045626-6348ff80-c1f8-11ea-9142-63735224f289.png" width="350" heigth="150">

5. Choose graphical layout. Users are allowed to choose between 5 different layouts.

- Circular layout.
- Force-directed layout (Fruchterman Reingold algorithm).
- Sphere layout.
- DH layout (Davidson-Harel algorithm). -only for network type whole-.
- Bipartite layout. -only for network type crosstalk-.

<img src="https://user-images.githubusercontent.com/49268378/87045880-b15e0300-c1f8-11ea-9ea2-10fd71ac9d48.png" width="150" height="125">

6. Enrichment Analysis

By default TALKIEN computes descriptive, topological network statistics and renders network plots. However, users are able to perform a functional analysis based on Reactome Pathway Database by selecting “Enrichment” tab, on top of the page. This may take between 10 seconds to 1 minute depending on the complexity of the lists.

<img src="https://user-images.githubusercontent.com/49268378/87046520-9b047700-c1f9-11ea-8c70-21476dce3b6a.png" width="150" height="50">

All results could be downloaded by clicking donwload buttons on the bottom of the tabs. Network plots could be downloaded first as html objects. Additionally, there is a black button to export the image to a png format once the downloaded html has been opened. 

<img src="https://user-images.githubusercontent.com/49268378/87046925-18c88280-c1fa-11ea-9bbd-f04d84545737.png" width="400" height="350">

MAIN TAB: PLOTS

Network will be displayed and able to be downloaded as an interactive html object. The size of the node is proportional to its degree and the width of the edge is proportional to its interaction score. In bipartite mode. There is the option of changing network orientation (UP-DOWN or LEFT-RIGHT) and changing the groups (bipartite network by lists or by location)

<img src="https://user-images.githubusercontent.com/49268378/87046571-a9529300-c1f9-11ea-96f5-be59aed316a9.png" width="250" height="150">

In addition, users can select one node by a selection list at top-left of the plot

<img src="https://user-images.githubusercontent.com/49268378/87046600-b1123780-c1f9-11ea-815d-c0927d8410b9.png" width="150" heigth="250">

Or by clicking on any node. If the selection is done by the latter option, a table with the node’s belonging pathways will be displayed. It shows pathway ID, pathway name and a link to see more in-depth pathway properties on reactome website. 

<img src="https://user-images.githubusercontent.com/49268378/87046632-b8d1dc00-c1f9-11ea-92f3-075524af383d.png" width="850" heigth="750">


NETWORK PARAMETERS TAB

For the plotted network, the network descriptive parameters are displayed: number of nodes and edges, diameter, shortest_Pathway, density, average neighbors, clustering_coefficient and centralization.

For all nodes five different centralization methods are displayed and can be downloaded. For each pair of connected nodes, interaction score, location, list and other annotation IDs are displayed and can be downloaded.
ENRICHMENT ANALYSIS TAB

In this tab are displayed the enrichment analysis results in two different ways.

- Table summarizing enriched pathways, proteins involved, p-values, enrichment scores, ect.
- Dotplot with up to 15 top enriched pathways.

