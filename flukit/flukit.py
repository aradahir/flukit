#!/usr/bin/env python3

import sys
import typer
import pandas

from pathlib import Path
from rich import print
from Bio import SeqIO
from flukit.utils.variants import set_gene
from flukit.utils.run import call_variants, call_clades

app = typer.Typer(
    help = "flukit - the influenza surveillance toolkit... kinda",
    add_completion=False)

@app.command(no_args_is_help=True)
def main(
    sequences: Path = typer.Option(
        ...,
        "-s",
        "--sequences",
        help="Path to sequences."), 
    lineage: str = typer.Option(
        ...,
        "-l",
        "--lineage",
        help="lineage options are: h1n1, h3n2, vic."),
    output: Path = typer.Option(
        ...,
        "-o",
        "--output",
        help="output path"

    ),
    batchNumber: str = typer.Option(
        None,
        "-b",
        "--batchNumber",
        help="prefix used for output files, optional.")
     ):

    # input checks
    if lineage not in ['h1n1', 'h3n2', 'vic']:
        raise typer.BadParameter(f"Linease '{lineage}' is not valid. You must choose from: \n h1n1\n h3n2\n vic")
    if not Path(sequences).resolve():
        raise typer.BadParameter(f"The path is not correct, please check: {sequences}")

    # outputs
    if batchNumber is not None:
        results_out = Path(output) / (batchNumber + "_results.csv")
        clades_out = Path(output) / (batchNumber + "_clades.txt")
    else:
        results_out = Path(output) / "results.csv"
        clades_out = Path(output) / "clades.txt"

    # parse input sequences
    try:
        input_sequences = SeqIO.to_dict(SeqIO.parse(sequences, "fasta"))
        input_sequences = set_gene(input_sequences)
    except FileNotFoundError:
        raise typer.BadParameter(f"Check input: {sequences} \nFile does not exist.")
    except ValueError as error:
        raise typer.BadParameter(f'''Check input: {sequences} \nError: {error}''')
    except:
        raise typer.BadParameter(f"[bold yellow] Error reading in fasta file. Check input: {sequences} \nError: {sys.exc_info()[0]}")

    # call variants
    variants, ha_records = call_variants(input_sequences, lineage)
    
    # call clades
    clades = call_clades(ha_records, lineage)
    
    # combine dataframes
    results = pandas.merge(
        variants, clades, 
        how="left", left_index=True, 
        right_on="seqName", 
        suffixes=(False, False)
        )
    results.insert(0, 'seqName', results.pop('seqName')) # reorder 
    
    # write results
    results.to_csv(results_out, sep=',', index=False)
    print("[bold green]All done![/bold green]")
