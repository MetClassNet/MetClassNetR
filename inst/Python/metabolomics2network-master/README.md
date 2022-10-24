# Metabolites2Network 1

The aim of the Python Metabolites2Network package is to provide methods and tools for a flexible way of matching metabolites identified using metabolomics or lipidomics to networks.

This flexible matching of identifiers is based on ontology. 
In particular, we use [ChEBI ontology](https://www.ebi.ac.uk/chebi/) to allow matching of precisely identified molecules (e.g. lipid species) to more generic descriptions of molecules that can be found in metabolic networks (e.g. lipid classes).

The package works for metabolomics data but is particularly well suited for lipids.

The program takes as input the metabolites/lipids annotated and the metabolites belonging to a network.
Both have to use ChEBI as identifiers.
Then the program returns closest matches using the ChEBI ontology Directed Acyclic Graph.

Metabolites2Network code is provided by INRA MetExplore group and available in [MetExplore web server](http://www.metexplore.fr/).
This project was developed within [MetaboHub: France metabolomics and fluxomics infrastructure](https://www.metabohub.fr/] and INRA) and [TOXALIM](https://www6.toulouse.inra.fr/toxalim_eng/) laboratory.

## Mapping result examples

These examples are computed on  Recon 2.2 metabolic network (BioSource 4311 in MetExplore, coming from [Swainston et al. 2016](https://www.ncbi.nlm.nih.gov/pubmed/27358602) )


| ChEBI id in Dataset        | Dataset molecule name         | ChEBI id in Recon2.2 | matching molecules in Recon2.2  | ontology connection  | mapping distance  |
| ------------- |:-------------:| -----:|:-------------:|:-------------:|:-------------:|
| CHEBI:15756      |  hexadecanoic acid | 	CHEBI:7896 | Palmitate (hexadecanoate), M_hdca | hexadecanoic acid (CHEBI:15756) **is conjugate acid of** hexadecanoate (CHEBI:7896) | -0.1 |
| CHEBI:90488      |  phosphatidylinositol (18:1/20:4) | 	CHEBI:57880 | 1-phosphatidyl-1D-myo-inositol(1−), M_pail_hs | phosphatidylinositol (18:1/20:4) (CHEBI:90488) **is a** 1-phosphatidyl-1D-myo-inositol (CHEBI:16749) <br/><br/> 1-phosphatidyl-1D-myo-inositol (CHEBI:16749) **is conjugate acid of** 1-phosphatidyl-1D-myo-inositol(1−) (CHEBI:57880) | 1.1 |
| CHEBI:36023      |  vaccenic acid | 	CHEBI:30828 | trans-vaccenate, M_vacc | trans-vaccenate(1−) (CHEBI:30828) **is conjugate base of** trans-vaccenic acid (CHEBI:28727) </br></br> trans-vaccenic acid (CHEBI:28727) **is a** vaccenic acid (CHEBI:36023)| -1.1 |

Here are some explanations of these results:

- First line shows a charge difference between the metabolite in the dataset and the one present in the network.
- Second line shows an example where the metabolite in the dataset corresponds to a species while the one in the network corresponds to the class. This results in a positive distance (1.1).
- Third line shows an example where the metabolite in the dataset is a class while the metabolite definition is more precise in the network. This results in a negative distance (-1.1).


## Setting up python and installing required libraries

### 1. Install python
- Download and python from [https://www.python.org/downloads/](https://www.python.org/downloads/) make sure you get a 3.X version and not a 2.X, the scripts were written using python 3.4.6
- During python installation check the box : "Add python to PATH" so that you can easily call the scripts later
- During the installation also install pip that will be used to download and install the required libraries in the next steps [https://pypi.org/project/pip/](https://pypi.org/project/pip/)

### 2. Install required libraries:
- In your console, set the current directory to the directory containing all the python files
- In command line run :

```pip install -r requirements.txt -U```
	

## Overview of the command and parameters

### Command line

Command line works as follows: 

```python3 ./metabolomics2network.py file_type metabolomics_path network_metabolites_path output_path json_conf_file_path```

### Arguments

**file_type**

First parameter specify the type of files used as input for metabolomics data and for output file.
Values has to be ```json``` or ```tsv```.

**metabolomics_path:** Path to the metabolomics (including lipidomics) dataset

This file contains the metabolites that will be matched to the network. Here is a small example:

```
[{"name":"M1","undefined":"","chebi":"17408"}]
```

**network_metabolites_path:** Path to the network metabolites json file

The file contains all the metabolites of a given metabolic network formatted in json. Here is a small example:

```
[{"id":"7230990","name":"6-hydroxypaclitaxel","dbIdentifier":"M_htaxol_b","chemicalFormula":"C47H51NO15","idin":[{"extDBName":"chebi","extID":"CHEBI:63859","origin":"SBML File","score":"1"},{"extDBName":"inchi","extID":"InChI=1S\/C47H51NO15\/c1-24-30(61-43(57)33(51)32(27-16-10-7-11-17-27)48-41(55)28-18-12-8-13-19-28)22-47(58)40(62-42(56)29-20-14-9-15-21-29)36-45(6,38(54)35(60-25(2)49)31(24)44(47,4)5)37(53)34(52)39-46(36,23-59-39)63-26(3)50\/h7-21,30,32-37,39-40,51-53,58H,22-23H2,1-6H3,(H,48,55)\/t30-,32-,33+,34-,35+,36-,37-,39+,40-,45-,46+,47+\/m0\/s1","origin":"SBML File","score":"1"}]},
{"id":"7230991","name":"6-hydroxypaclitaxel","dbIdentifier":"M_htaxol_c","chemicalFormula":"C47H51NO15","idin":[{"extDBName":"chebi","extID":"CHEBI:63859","origin":"SBML File","score":"1"},{"extDBName":"inchi","extID":"InChI=1S\/C47H51NO15\/c1-24-30(61-43(57)33(51)32(27-16-10-7-11-17-27)48-41(55)28-18-12-8-13-19-28)22-47(58)40(62-42(56)29-20-14-9-15-21-29)36-45(6,38(54)35(60-25(2)49)31(24)44(47,4)5)37(53)34(52)39-46(36,23-59-39)63-26(3)50\/h7-21,30,32-37,39-40,51-53,58H,22-23H2,1-6H3,(H,48,55)\/t30-,32-,33+,34-,35+,36-,37-,39+,40-,45-,46+,47+\/m0\/s1","origin":"SBML File","score":"1"}]},]
```

The file provided as an example in the repository was produced using MetExplore web server (wwww.metexplore.fr).
This file corresponds to the metabolites of Recon 2.2 metabolic network (BioSource 4311 in MetExplore, coming from [Swainston et al. 2016](https://www.ncbi.nlm.nih.gov/pubmed/27358602) )
		
If you need a specific file please contact: contact-metexplore@inra.fr 
			``
**output_path:** Path where output files will be writen

Depending on the ```file_type``` option, the file will be returned in txt (tsv) or json

**conf_file_path:** Path to the file containing aliases to ensure correspondance in identifiers labels.
Here is an example:

```
name alias: name
chebi alias: chebi
formula alias: formula
kegg alias: kegg
pubchem alias: pubchem
hmdb alias: hmdb
inchi alias: inchi
smiles alias: smiles
lipidmaps alias: lipidmap
swisslipids alias: swisslipids
```

**mapping_types:** 
- 1 = exact multimapping
- 2 = chebi class mapping


*Remarks:* 
- *To use multiple mapping at once use commas, for example 1,2 for exact mapping and ChEBI mapping*
- *Output consider results in this order of priority: Exact multimapping > ChEBI class mapping*
	
## Mapping on ChEBI classes using json files (metabolites of network and lipidomic results)	

```python3 ./metabolomics2network.py json ./example_data/data_15756.json ./example_data/metabolites_4311_DB.json ./example_data/data_out.json ./conf.txt 1,2```

Use the ```json``` parameter.

### Metabolomics/lipidomics json data file

This is a basic example file for data:

```
[{"name":"M1","undefined":"","chebi":"15756"}]
```

### Output in json format

Here is an example of output:

```
[
    {
        "name": "M1",
        "undefined": "",
        "chebi": "15756",
        "mapped": [
            {
                "path": [
                    "15756"
                ],
                "datasetname": "M1",
                "idsql": [
                    "7236116"
                ],
                "networkid": "M_hdca_x",
                "mapping_type": "chebi class mapping",
                "distance": "-0.1",
                "identifiers": {
                    "chebi": [
                        "7896",
                        "7896"
                    ]
                }
            },
            {
                "path": [
                    "15756"
                ],
                "datasetname": "M1",
                "idsql": [
                    "7236408"
                ],
                "networkid": "M_hdca_e",
                "mapping_type": "chebi class mapping",
                "distance": "-0.1",
                "identifiers": {
                    "chebi": [
                        "7896",
                        "7896"
                    ]
                }
            },
            {
                "path": [
                    "15756"
                ],
                "datasetname": "M1",
                "idsql": [
                    "7236411"
                ],
                "networkid": "M_hdca_c",
                "mapping_type": "chebi class mapping",
                "distance": "-0.1",
                "identifiers": {
                    "chebi": [
                        "7896",
                        "7896"
                    ]
                }
            },
            {
                "path": [
                    "15756"
                ],
                "datasetname": "M1",
                "idsql": [
                    "7236416"
                ],
                "networkid": "M_hdca_l",
                "mapping_type": "chebi class mapping",
                "distance": "-0.1",
                "identifiers": {
                    "chebi": [
                        "7896",
                        "7896"
                    ]
                }
            },
            {
                "path": [
                    "15756"
                ],
                "datasetname": "M1",
                "idsql": [
                    "7236417"
                ],
                "networkid": "M_hdca_r",
                "mapping_type": "chebi class mapping",
                "distance": "-0.1",
                "identifiers": {
                    "chebi": [
                        "7896",
                        "7896"
                    ]
                }
            },
            {
                "path": [
                    "15756"
                ],
                "datasetname": "M1",
                "idsql": [
                    "7236435"
                ],
                "networkid": "M_hdca_b",
                "mapping_type": "chebi class mapping",
                "distance": "-0.1",
                "identifiers": {
                    "chebi": [
                        "7896",
                        "7896"
                    ]
                }
            }
        ]
    }
]
```

	
## Mapping on ChEBI classes using tsv files (metabolites of network and lipidomic results)

1. open the conf_tsv.txt file in the main directory to set the names of the columns in the metabolites file and the lipid file
There is no need for every column to be present in your files
In case some metabolites can have multiple occurences of the same type of identifiers, you can also set the separator character (by default "|")

2. run the mapping script using command line:

```python3 ./metabolomics2network.py tsv ./example_data/data_15756.txt ./example_data/metabolites_4311_DB.json ./example_data/data_out.txt ./conf.txt 1,2```

The ```tsv``` option has to be set.

Data has to be formatted following order specified in the conf.txt file and **has to start with an header line**.
Here is an example:

```
name	chebi
hexadecanoic acid	15756
```

Out file is returned in tsv format. For example:

```
metabolite name	mapped on id	mapping types	distance	chebi	inchi	smiles	hmdb	kegg	pubchem	lipidmaps	swisslipids
hexadecanoic acid	M_hdca_x;M_hdca_e;M_hdca_c;M_hdca_l;M_hdca_r;M_hdca_b;	chebi class mapping;chebi class mapping;chebi class mapping;chebi class mapping;chebi class mapping;chebi class mapping;	-0.1;-0.1;-0.1;-0.1;-0.1;-0.1;	['7896', '7896'];['7896', '7896'];['7896', '7896'];['7896', '7896'];['7896', '7896'];['7896', '7896'];	N/A;N/A;N/A;N/A;N/A;N/A;	N/A;N/A;N/A;N/A;N/A;N/A;	N/A;N/A;N/A;N/A;N/A;N/A;	N/A;N/A;N/A;N/A;N/A;N/A;	N/A;N/A;N/A;N/A;N/A;N/A;	N/A;N/A;N/A;N/A;N/A;N/A;	N/A;N/A;N/A;N/A;N/A;N/A;
```
						   
## Integration testing

### Testing json based mapping

```python3 ./tests/test_integration_mapping_json.py```

## Release Notes

**v1.0.2** updated the package to fit release 1.0.9 of libChEBIpy

## Built With

* [LibChEBIpy](https://github.com/libChEBI/libChEBIpy) LibChEBIpy (see [Swainston et al. 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4772646/) for more details)

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Arthur Moreau** - *Initial code*
* **Florence Vinson** - *Debugging and feature add-ons*
* **Fabien Jourdan** - *Debugging and feature add-ons*

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Test, controls and use cases were provided by Nathalie Poupin
* This project is part of MetaboHub French National Metabolomics and Fluxomics Infrastructure
* This project is also an outcome of Metatoul Toulouse metabolomics and fluxomics facility
* Project was funded by INRA Halomics project (PI Dr. Nicolas Cabaton)

## Partners
![MetExplore](/images/MetExplore.png)

![MetaboHub](/images/MetaboHub.png)

![Metatoul](/images/Metatoul.jpg)