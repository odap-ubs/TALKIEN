<img src="https://user-images.githubusercontent.com/49268378/87043206-ff710780-c1f4-11ea-81c5-36523124ddf0.jpg" width="400" height="100"> 
TALKIEN - crossTALK bIpartite Network -

You can plot customized networks, obtain network parameters and enrichment analysis.
INPUT

User is required to input two list of genes/proteins. Each of them belongs to one of the two entities the user in interested in. To decipher the molecular crosstalk between them, TALKIEN only uses secreted and membrane receptors proteins. The hypothesis behind is that cellular communication could be partially explained by physical interactions between receptor proteins and their specific ligands. Specifically in paracrine signalling, receptors are activated by ligands proteins secreted by another cell. However, TALKIEN also allow users to visualize pathways downstream receptors in the network (see below): 

<img src="https://user-images.githubusercontent.com/49268378/87044332-ac984f80-c1f6-11ea-9c6b-1c2189a5a01f.jpeg" width="350" height="150">

These lists may derive from a variety of experiments and could be measured with diverse transcriptomics (such as RNA-seq, microarray or scRNA-seq) and proteomics approaches (mass spectrometry or protein chips among others). For example, co-culture experiments of cancer cells with/without stromal cells, or cells growing under the stimulation of supernatant media collected from immune cells cultures. Lists from the example came from the following transcriptomic experiment: 

<img src="https://user-images.githubusercontent.com/49268378/87044646-216b8980-c1f7-11ea-840a-53570c701504.png" witdh="350" height="350">

PLOT NETWORKS

The basic elements of the networks are nodes (Proteins) and edges (interactions between them) that connect nodes within the network. Connections are defined according to a score based on certain evidence of interaction.
User options.

1. Upload lists Input List #1 and Input List #2. Upload files with header and one column lenght. Otherwise column #1 will be used as reference. Both files must have the same annotation: entrez ID, ensembl ID, gene symbol or uniprot ID are allowed.
<img src="https://user-images.githubusercontent.com/49268378/87044850-68f21580-c1f7-11ea-901f-55e3175934d5.png" width="350" height="150">

There is a preloaded example with two lists. By selecting “Load example data” option, an example network will be displayed. 
<img src="https://user-images.githubusercontent.com/49268378/87045008-9fc82b80-c1f7-11ea-9be9-db4b6ec6a971.png" width="200" height="50">


2. Choose annotation type.

- Gene Symbol
- Entrez
- Ensembl
- Uniprot

<img src="https://user-images.githubusercontent.com/49268378/87045210-e74eb780-c1f7-11ea-8939-699b63014aee.png" width="250" height="150">

3. Choose network type. Users are allowed to compute interactions between lists (the so-called crosstalk) or interactions between and within lists.

- whole –> get interactions from list 1 to list1, from list2 to list2, from list1 to list2, and viceversa
- crosstalk –> get interactions only from list1 to list2 and viceversa

<img src="https://user-images.githubusercontent.com/49268378/87045423-27ae3580-c1f8-11ea-8030-d618b0b81860.png" width="150" height="150">

4. Select score threshold. That is the minimum interaction value to link two proteins. Only protein-protein interaction with scores higher than chosen will be shown in the network. Score threshold is based on STRING’s combined score. 

<img src="https://user-images.githubusercontent.com/49268378/87045626-6348ff80-c1f8-11ea-9142-63735224f289.png" width="150" heigth="150">

5. Choose graphical layout. Users are allowed to choose between 5 different layouts.

- Circular layout.
- Force-directed layout (Fruchterman Reingold algorithm).
- Sphere layout.
- DH layout (Davidson-Harel algorithm). -only for network type whole-.
- Bipartite layout. -only for network type crosstalk-.

