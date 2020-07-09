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

    Upload lists Input List #1 and Input List #2. Upload files with header and one column lenght. Otherwise column #1 will be used as reference. Both files must have the same annotation: entrez ID, ensembl ID, gene symbol or uniprot ID are allowed.


