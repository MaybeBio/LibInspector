import typer
from typing import Optional
from .lib_inspector import inspect_library

# optimize: Instead of defining a Typer app and subcommands,
# define a core function that serves directly as the CLI entry point

def main_action(
    library: str = typer.Argument(..., help="Name of the library to inspect (e.g. 'numpy', 'Bio.PDB')."),
    output: Optional[str] = typer.Option(None, "--output", "-o", help="Path to save the Markdown report. If not provided, prints to stdout."),
    private: bool = typer.Option(False, "--private", help="Include private members (starting with '_')."),
    imported: bool = typer.Option(False, "--imported", help="Include members imported from other modules."),
    limit_api: int = typer.Option(20, "--limit-api", help="Limit for API recommendations."),
    limit_snippets: int = typer.Option(20, "--limit-snippets", help="Limit for code snippets."),
    limit_args: int = typer.Option(3, "--limit-args", help="Threshold for multi-line arguments in snippets."),
    limit_ext: int = typer.Option(20, "--limit-ext", help="Limit for external libraries."),
    limit_pr: int = typer.Option(20, "--limit-pr", help="Limit for PageRank metrics."),
    limit_dep: int = typer.Option(100, "--limit-dep", help="Limit for dependency graph nodes."),
    limit_inh: int = typer.Option(100, "--limit-inh", help="Limit for inheritance graph classes."),
    limit_guide: int = typer.Option(20, "--limit-guide", help="Limit for extraction guide functions.")
):
    """
    Inspect a Python library/module and generate Markdown/Html documentation.
    
    \b
    Example:
    \b
    # Inspect a CLI tool library (e.g. metapredict-predict-disorder --help)
    lib-inspector metapredict-predict-disorder --private -o Metapredict_CLI_Docs.md
    \b
    # Inspect a Standard library module (e.g. import numpy)
    lib-inspector numpy --imported -o Numpy_Docs.md
    \b
    # Inspect YOUR local library
    lib-inspector ./my_project/train.py -o ./reports/train.md

    \b
    Notes:
    - 1, --private means to include private members, while by default they are excluded if not provided.
    - 2, --imported means to include members imported from other modules, while by default they are excluded if not provided.
    """
    inspect_library(library_name=library, 
                    output_path=output, 
                    include_private=private, 
                    include_imported=imported,
                    limit_api_recommendations=limit_api,
                    limit_code_snippets=limit_snippets,
                    limit_snippet_args=limit_args,
                    limit_external_libs=limit_ext,
                    limit_pagerank=limit_pr,
                    limit_dependency_graph=limit_dep,
                    limit_inheritance_graph=limit_inh,
                    limit_extraction_guide=limit_guide)

def main():
    # Use typer.run to make the function the main command
    typer.run(main_action)

if __name__ == "__main__":
    main()