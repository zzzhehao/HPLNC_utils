## HPLNC

Utilities and handy functions I wrote for my Master Thesis.

## Installation

Install from Github:

```r
devtools::install_github("https://github.com/zzzhehao/HPLNC_utils.git")
```

## Function

Some functions are for general purpose and could easily be used in other projects, the others are more or less project-based, adapted to the file structure I am using. They are mostly nice wrapper functions for handling repetitive work in my daily data management. Some could be easily modified to a general purpose function (but I am not going to do it soon because it's extra work for me), some are very specific in file handling so that it makes no sense to modify them into general purpose functions, but maybe you can get some inspiration from them.

| **Function**             | **Purpose**                        | **Description**                                                                              | **Project-specific**    |
| ------------------------ | ---------------------------------- | -------------------------------------------------------------------------------------------- | ----------------------- |
| `DBchange_sign`          | Database Utility                   | Generate Signature for Database Change                                                       | Yes                     |
| `DBexecute`              | Database Utility                   | Execute YAML Change Request                                                                  | Yes                     |
| `DBlastRequest`          | Database Utility                   | Return Signature of the Last Request                                                         | Yes                     |
| `DBmanual_pull`          | Database Utility                   | Pull Raw Table                                                                               | Yes                     |
| `DBmanual_write`         | Database Utility                   | Write Table                                                                                  | Yes                     |
| `DBpullTable`            | Database Utility                   | Pull Table                                                                                   | Yes                     |
| `DBsnapshot`             | Database Utility                   | Generate Snapshot from All Tables                                                            | Yes                     |
| `create_log`             | File Utility                       | Initialize daily log                                                                         | Yes                     |
| `date2assets`            | File Utility                       | Create assets folder for the day                                                             | Yes                     |
| `dms_to_decimal`         | GIS                                | Transform DMS coordinates into DD format.                                                    | No                      |
| `dec_format`             | GIS                                | Format DD coordinates into pure numerical object.                                            | No                      |
| `voucher.inquery`        | Metadata Utility                   | Query specimen info with voucher                                                             | Yes                     |
| `ID.inquery`             | Metadata Utility                   | Query specimen info with DZMB2HH ID                                                          | Yes                     |
| `load_morphocheck`       | Metadata Utility                   | Load latest morphocheck result                                                               | Yes                     |
| `summary_sp`             | Metadata Utility                   | Generate Species Summary                                                                     | Yes                     |
| `update_INSD_metadata`   | Metadata Utility, Database Utility | Wrapper Function to Read All XML File and Update Information to Databse                      | Yes                     |
| `update_sequence_map`    | Metadata Utility, Database Utility | Create a Summary Table Mapping Sequence Label/Accession Number with Vouchered Animals.       | Minimal, easy to modify |
| `bathy.basemap`          | Plotting                           | Wrapper function for generating bathymetric + land elevation basemap in ggplot for your map. | No                      |
| `geom_checkerboard`      | Plotting                           | Create checkerboard around the map (only for mercator projection).                           | No                      |
| `midwayTreeViz`          | Plotting                           | Snapshot the MrBayes Result During the Analysis                                              | Minimal, easy to modify |
| `visualize_tree`         | Plotting                           | Visualize Consensus Tree Produced From MrBayes Analysis                                      | Minimal, easy to modify |
| `fetch.land.elv`         | Plotting, GIS                      | Fetch land elevation data for plotting purpose.                                              | No                      |
| `fetch.NOAA.bathy`       | Plotting, GIS                      | Fetch NOAA bathymetric data for plotting purpose.                                            | No                      |
| `create_summary_sp`      | Plotting, Metadata Utility         | Create Species Summary Report in Image and PDF                                               | Yes                     |
| `gblocks`                | Program Wrapper                    | Run Gblocks to Trim the Alignment                                                            | No                      |
| `render_species_summary` | Program Wrapper                    | Wrapper for Quarto Render Species Summary                                                    | Yes                     |
| `clearRcache`            | R Utility                          | Clear temporal file in R session.                                                            | No                      |
| `clean_label`            | Sequence Utility                   | Clean Sequence Labels in Alignment                                                           | Yes                     |
| `extract_INSD_metadata`  | Sequence Utility                   | Extract INSD Metadata from xml File Downloaded From NCBI Genebank                            | No                      |
| `nexus2fasta`            | Sequence Utility                   | Convert NEXUS alignment to FASTA alignment                                                   | No                      |

All functions have documentation that can be called with `help()` or `?` in R. Also check the PDF reference manual for details.

Thanks Gemini for keeping some function tidy and robust when I was exhausted. 
