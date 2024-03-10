# Advanced-flowcytoscript
An optimized high dimensional workflow for flow cytometry data designed for people comfortable with some R.

Code is available and free for academic users. Commercial users should contact Adrian Liston to discuss licensing options.

This version of the flowcytoscript has similar features as the Simplified-flowcytoscript, but without the interactive workflow. All aspects of the analysis are customizable, as are all the graphical output parameters. It's easier to repeat sections and re-run the analysis with modifications.

To select channels for analysis, use the get_channels script.

Initial installation of packages can be automated with the installation script.

Flowcytoscript will automatically detect the cytometer used for collection of the data and apply an appropriate biexponential transformation. This will reduce errors and artefacts due to incorrect scaling (these can be big problems).

If you wish to set the scaling explicitly, we recommend using the CSV-based method. This allows users to set the data scales (transformations) in FlowJo or any other standard flow cytometry data analysis software. Scaling of the data is critical for optimal visualization and clustering. To understand how to set scales in FlowJo, see the instruction file "Setting axis in FlowJo for Aurora data.pptx". Use the channel values format. To export your data in CSV format, preserving the transformations from FlowJo, see the instructions in "Exporting data in csv format.PNG". More detail [here](see https://docs.flowjo.com/flowjo/graphs-and-gating/gw-transform-overview/)
