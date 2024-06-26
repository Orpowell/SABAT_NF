#!/usr/bin/env python3
import logging
import sys
import click
from blast_2_bed import blast2bed
from assemble_gene import assemble_exons
from assemble_gene import assemble_locus

logging.basicConfig(
    stream=sys.stderr,
    format="%(asctime)s - %(message)s",
    datefmt="%d-%b-%y %H:%M:%S",
    level=logging.INFO,
)


@click.group(help="S.A.B.A.T: Semi-Automatic BLAST Annotation Toolkit")
@click.version_option("-v", "--version", message="SABAT 0.3.0")
def cli() -> None:
    pass


cli.add_command(blast2bed)
cli.add_command(assemble_exons)
cli.add_command(assemble_locus)

if __name__ == "__main__":
    cli()
