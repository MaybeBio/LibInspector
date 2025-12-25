import inspect # Inspect live objects and class internals
import importlib # Dynamically import modules
import sys
import os
import pkgutil # Walk through packages/modules (used to find submodules)
import ast # Parse Python code into Abstract Syntax Trees (AST) for code analysis/refactoring
import re
from collections import Counter, defaultdict, deque
from typing import Any, List, Dict, Optional, Tuple, Set
import json
import html # Used for Mermaid MD to HTML conversion, escaping HTML characters
from pathlib import Path # For file path handling


# Import metadata for CLI/Package resolution
try:
    from importlib import metadata as importlib_metadata
except ImportError:
    import importlib_metadata # Fallback for older Python versions


# Attempt to import networkx for advanced network analysis
try:
    import networkx as nx
    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False


# --- Helper: ID Sanitizer ---
def sanitize_id(name: str) -> str:
    """
    Description
    ----------
    Convert an arbitrary string into a valid Mermaid node ID. Mermaid node IDs
    should only contain letters, digits and underscores.

    Args
    -----
    name : str
        The input string to sanitize into a Mermaid-safe identifier.

    Returns
    --------
    str
        A sanitized identifier containing only [A-Za-z0-9_] and guaranteed not
        to start with a digit (an underscore is prepended if necessary).

    Notes
    -------
    - 1, Non-alphanumeric characters (including dots, spaces and symbols) are
      replaced with underscores.
    - 2, If the resulting identifier starts with a digit, a leading underscore is
      added to ensure it is a valid ID for Mermaid.
    - 3, Empty input returns an empty string.
    """
    # Replace dots, Spaces, and special symbols with underscores
    clean = re.sub(r'[^a-zA-Z0-9_]', '_', name)
    # Avoid starting with a digit
    if clean and clean[0].isdigit():
        clean = "_" + clean
    return clean


# --- Helper: Target Resolver (Intelligent Input Handling) ---
class ResolvedTarget:
    """
    Description
    ----------
    Data class to hold the result of target resolution.
    """
    def __init__(self, import_name, target_type, original_input, file_path=None, cli_info=None, dist_name=None):
        self.import_name = import_name # The string to pass to importlib.import_module
        self.target_type = target_type # 'module', 'file', 'cli_command', 'distribution'
        self.original_input = original_input
        self.file_path = file_path # If it was a local file
        self.cli_info = cli_info # {'command': str, 'func_ref': str} if resolved from CLI
        self.dist_name = dist_name # The PyPI distribution name (e.g., 'scikit-learn')

class TargetResolver:
    """
    Description
    ----------
    Resolves ambiguous user inputs (CLI commands, file paths, PyPI names) into importable module names.
    """
    @staticmethod
    def resolve(target: str) -> ResolvedTarget:
        # 1. Check if it is a local file path
        path_obj = Path(target)
        if path_obj.exists() and (path_obj.is_file() or path_obj.is_dir()):
            abs_path = path_obj.resolve()
            if abs_path.is_file() and abs_path.suffix == '.py':
                # It's a python script
                module_name = abs_path.stem
                return ResolvedTarget(module_name, 'file', target, file_path=str(abs_path))
            elif abs_path.is_dir():
                # It's a package directory (containing __init__.py)
                if (abs_path / "__init__.py").exists():
                    module_name = abs_path.name
                    return ResolvedTarget(module_name, 'file', target, file_path=str(abs_path.parent))
        
        # 2. Check if it is a registered CLI command (console_scripts)
        # Example: 'metapredict-predict-disorder' -> maps to 'metapredict' module
        try:
            # Get all console_scripts
            eps = importlib_metadata.entry_points(group='console_scripts')
            # Filter matches
            matches = [ep for ep in eps if ep.name == target]
            if matches:
                ep = matches[0]
                # ep.module is the module containing the entry point (e.g., metapredict.scripts.predict_disorder)
                # We want the top-level package (e.g., metapredict)
                top_package = ep.module.split('.')[0]
                return ResolvedTarget(
                    import_name=top_package,
                    target_type='cli_command',
                    original_input=target,
                    cli_info={'command': target, 'func_ref': ep.value, 'module': ep.module},
                    dist_name=ep.dist.name if hasattr(ep, 'dist') and ep.dist else None
                )
        except Exception: pass

        # 3. Check if it is a Distribution Package name (PyPI name) that differs from import name
        # Example: 'scikit-learn' -> imports as 'sklearn'
        try:
            dist = importlib_metadata.distribution(target)
            # Try to find top_level.txt which lists the importable modules
            top_level = dist.read_text('top_level.txt')
            if top_level:
                imports = [x for x in top_level.splitlines() if x]
                if imports:
                    # Use the first top level module found
                    return ResolvedTarget(imports[0], 'distribution', target, dist_name=target)
        except importlib_metadata.PackageNotFoundError:
            pass

        # 4. Default: Assume it is a direct module name
        return ResolvedTarget(target, 'module', target)


# --- Helper: 1. Single-function Logic Analysis (Micro) ---
class LogicNode:
    """
    Description
    ----------
    Models a single element in a data/logic flow diagram extracted from a
    function's AST. Nodes denote inputs, processing steps, or outputs and
    are used to assemble Mermaid flowcharts showing data production and use.

    Attributes
    ----------
    id : str
        Unique identifier for the node (used when rendering the diagram).
    label : str
        Human-readable label displayed inside the node.
    node_type : str
        One of "input", "process", or "output" — controls node shape/style.
    edges_in : List[Tuple[str, str]]
        Incoming edges as (source_node_id, variable_name) pairs indicating
        which node produced each input variable.
    """
    def __init__(self, id, label, node_type="process"):
        self.id = id
        self.label = label
        self.node_type = node_type # input, process, output
        self.edges_in = [] # List of (source_id, var_name)

class AdvancedFlowVisitor(ast.NodeVisitor):
    """
    Description
    ----------
    Parses function source code to build a data flow graph.
    Tracks the chain of variable production (Definition) -> consumption (Usage).
    """
    def __init__(self):
        """
        Description
        ----------
        Initialize visitor state used to accumulate flow nodes and track current producers.

        Args
        -----
        None

        Returns
        --------
        None

        Notes
        -------
        - Initializes:
          - self.nodes: list of LogicNode instances discovered
          - self.current_producers: mapping var_name -> node_id for latest producer
          - self.counter: integer counter for generating unique node ids
        """
        self.nodes = []
        self.current_producers = {} # var_name -> node_id (Tracks which node currently produces each variable)
        
        self.counter = 0

    def _get_id(self):
        """
        Description
        ----------
        Generate a new unique internal node identifier.

        Args
        -----
        None

        Returns
        --------
        str
            A new unique node id string (e.g., "Node1").

        Notes
        -------
        - Increments an internal counter on each call.
        """
        self.counter += 1
        return f"Node{self.counter}"

    def _resolve_inputs(self, input_vars: List[str]) -> List[Tuple[str, str]]:
        """
        Description
        ----------
        Resolve which previously created nodes produced the given input variables.

        Args
        -----
        input_vars : List[str]
            Variable names referenced on the right-hand side of an expression.

        Returns
        --------
        List[Tuple[str, str]]
            List of (source_node_id, var_name) pairs for known producers.

        Notes
        -------
        - Only returns producers that were previously recorded in self.current_producers.
        """
        edges = []
        for var in input_vars:
            if var in self.current_producers:
                source_id = self.current_producers[var]
                edges.append((source_id, var))
        return edges

    def _extract_names(self, node) -> List[str]:
        """
        Description
        ----------
        Extract variable names referenced inside an AST node.

        Args
        -----
        node : ast.AST
            The AST node to inspect for Name and Attribute usage.

        Returns
        --------
        List[str]
            Unique list of variable names found (includes simple names and 'self.attr' for attributes).

        Notes
        -------
        - Only collects names used in load (read) context.
        - Captures simple `Name` nodes and `self.attr` attributes; other attributes are visited recursively.
        """
        names = []
        class NameCollector(ast.NodeVisitor):
            def visit_Name(self, n):
                if isinstance(n.ctx, ast.Load):
                    names.append(n.id)
            def visit_Attribute(self, n):
                # Attempt to capture self.xxx
                if isinstance(n.value, ast.Name) and n.value.id == 'self':
                    names.append(f"self.{n.attr}")
                self.generic_visit(n)
        
        if node:
            NameCollector().visit(node)
        return list(set(names)) # Remove duplicates

    def visit_FunctionDef(self, node):
        """
        Description
        ----------
        Handle a FunctionDef AST node: register input parameter node and traverse body.

        Args
        -----
        node : ast.FunctionDef
            The function definition AST node.

        Returns
        --------
        None

        Notes
        -------
        - Creates an "Input" LogicNode if the function has parameters and marks parameters as produced
          by that Input node so subsequent assignments can reference them as inputs.
        - Continues traversal into the function body to collect assignments, calls and return nodes.
        """
        args = []
        arg_labels = []
        
        # Extract parameters and type annotations
        all_args = node.args.args + node.args.kwonlyargs
        if node.args.vararg: all_args.append(node.args.vararg)
        if node.args.kwarg: all_args.append(node.args.kwarg)

        for arg in all_args:
            var_name = arg.arg
            args.append(var_name)
            
            # Attempt to get type annotation
            ann = ""
            if arg.annotation:
                try:
                    if hasattr(ast, 'unparse'):
                        ann = ": " + ast.unparse(arg.annotation)
                    else:
                        ann = ": " + str(arg.annotation)
                except: pass
            # Use HTML line breaks to avoid issues with \n in Mermaid
            arg_labels.append(f"• {var_name}{ann}")
            
        if args:
            node_id = "Input"
            # Use HTML tags to enhance display (bold title and line breaks)
            label = "<b>Input Data</b><br/>" + "<br/>".join(arg_labels)
            logic_node = LogicNode(node_id, label, node_type="input")
            self.nodes.append(logic_node)
            
            # Register these variables as produced by the Input node
            for arg in args:
                self.current_producers[arg] = node_id
                # Also register self.arg (simplified handling for common __init__ pattern)
                if 'self' in args:
                    self.current_producers[f"self.{arg}"] = node_id
        
        # Continue traversing the function body
        for item in node.body:
            self.visit(item)

    def visit_Assign(self, node):
        """
        Description
        ----------
        Handle simple assignment AST nodes by delegating to the assignment handler.

        Args
        -----
        node : ast.Assign
            The assignment AST node.

        Returns
        --------
        None
        """
        self._handle_assign(node, node.targets)

    def visit_AnnAssign(self, node):
        """
        Description
        ----------
        Handle annotated assignment (PEP 526) such as `x: int = value`.

        Args
        -----
        node : ast.AnnAssign
            Annotated assignment AST node.

        Returns
        --------
        None
        """
        # Handle annotated assignment: x: int = value
        if node.value:
            self._handle_assign(node, [node.target], annotation=node.annotation)

    def _handle_assign(self, node, targets, annotation=None):
        """
        Description
        ----------
        Core handler for assignment-like statements. Analyzes RHS inputs, determines operation label,
        creates a LogicNode representing the assignment/call/op, and updates producer mapping.

        Args
        -----
        node : ast.AST
            The original AST node for context (Assign or AnnAssign).
        targets : List[ast.AST]
            Left-hand side target nodes.
        annotation : Optional[ast.AST]
            Optional annotation for the target (for AnnAssign).

        Returns
        --------
        None

        Notes
        -------
        - Recognizes Calls, BinaryOps and Constants on the RHS to provide a richer label.
        - Supports simple target types: Name and self.Attribute.
        - Updates self.current_producers so later statements can resolve dependencies.
        """
        input_vars = self._extract_names(node.value)
        
        # --- 1. Intelligent label generation and noise filtering ---

        label = "Assign"
        node_type = "process" # Default to normal processing node
        
        # A. Function call (key focus)
        if isinstance(node.value, ast.Call):
            func_name = self._get_func_name(node.value)
            
            # Filter out simple type conversions that do not help understand logic
            if func_name in ['int', 'str', 'float', 'list', 'tuple', 'len', 'print', 'type']:
                return 
            
            label = f"<b>Call:</b> {func_name}"
            
            # [Highlight] Core transformation recognition: AI models, mathematical operations, I/O operations
            # These are the "Black Box" internals that beginners care most about
            core_keywords = [
                'torch', 'np', 'numpy', 'pd', 'pandas', # Libraries
                'model', 'net', 'layer', 'encoder', 'decoder', # Model objects
                'predict', 'forward', 'transform', 'encode', 'decode', 'process', # Core actions
                'read', 'load', 'save', 'write', 'get', 'post' # I/O
            ]
            if any(k in func_name.lower() for k in core_keywords):
                node_type = "core_process"

        # B. Binary operations
        elif isinstance(node.value, ast.BinOp):
            op = type(node.value.op).__name__
            label = f"<b>Op:</b> {op}"
        
        # C. Constant assignment (usually configuration, beginners don't need to care unless it's a default value)
        elif isinstance(node.value, ast.Constant):
             return # Directly ignore to reduce chart noise
        
        # --- 2. Output variable processing (with type annotation) ---
        outputs = []
        output_labels = []
        for target in targets:
            if isinstance(target, ast.Name):
                var_name = target.id
                outputs.append(var_name)
                
                ann_str = ""
                if annotation and hasattr(ast, 'unparse'):
                    try: ann_str = ": " + ast.unparse(annotation).replace('\n', '')
                    except: pass
                output_labels.append(f"{var_name}{ann_str}")
            elif isinstance(target, ast.Attribute):
                if isinstance(target.value, ast.Name) and target.value.id == 'self':
                    var_name = f"self.{target.attr}"
                    outputs.append(var_name)
                    output_labels.append(var_name)

        if outputs:
            node_id = self._get_id()
            # Use HTML line breaks to separate "action" and "result" for better display
            full_label = f"{label}<br/>⬇<br/>" + ", ".join(output_labels)
            
            logic_node = LogicNode(node_id, full_label, node_type=node_type)
            logic_node.edges_in = self._resolve_inputs(input_vars)
            
            self.nodes.append(logic_node)
            
            for out in outputs:
                self.current_producers[out] = node_id

    # --- Control flow visualization (Decision Logic) ---
    def visit_If(self, node):
        """
        Description
        ----------
        Visualize decision points. Beginners need to know if there's a fork in the logic.
        """
        test_code = "Condition"
        if hasattr(ast, 'unparse'):
            try: test_code = ast.unparse(node.test).replace('\n', ' ').replace('"', "'")
            except: pass
            
        node_id = self._get_id()
        label = f"<b>Decision</b><br/>If {test_code}?"
        logic_node = LogicNode(node_id, label, node_type="decision")
        
        input_vars = self._extract_names(node.test)
        logic_node.edges_in = self._resolve_inputs(input_vars)
        self.nodes.append(logic_node)
        
        # Continue traversing child nodes
        self.generic_visit(node)
    
    
    
    def visit_Expr(self, node):
        """
        Description
        ----------
        Handle expression statements, commonly used for standalone function calls with side-effects.

        Args
        -----
        node : ast.Expr
            The expression AST node.

        Returns
        --------
        None

        Notes
        -------
        - Creates a LogicNode for the call but does not update producers because such calls usually
          don't assign to variables.
        """
        # Handle standalone function calls (no assignment), e.g., print(), model.eval()
        if isinstance(node.value, ast.Call):
            input_vars = self._extract_names(node.value)
            func_name = self._get_func_name(node.value)
            
            node_id = self._get_id()
            logic_node = LogicNode(node_id, f"Call: {func_name}")
            logic_node.edges_in = self._resolve_inputs(input_vars)
            
            self.nodes.append(logic_node)
            # Such calls usually have side effects but no explicit return variables, so current_producers is not updated

    def visit_Return(self, node):
        """
        Description
        ----------
        Handle return statements by creating an output LogicNode that links to its input expressions.

        Args
        -----
        node : ast.Return
            The return statement AST node.

        Returns
        --------
        None

        Notes
        -------
        - Attempts to stringify the returned expression when possible for clearer labels.
        """
        input_vars = []
        ret_str = "None"
        if node.value:
            input_vars = self._extract_names(node.value)
            if hasattr(ast, 'unparse'):
                try: 
                    # Remove newline characters to prevent breaking Mermaid syntax
                    ret_str = ast.unparse(node.value).replace('\n', ' ')
                except: pass
            else:
                ret_str = "Expression"
        
        node_id = "Return"
        logic_node = LogicNode(node_id, f"Return\\n{ret_str}", node_type="output")
        logic_node.edges_in = self._resolve_inputs(input_vars)
        self.nodes.append(logic_node)

    def _get_func_name(self, node):
        """
        Description
        ----------
        Derive a human-readable function name for an ast.Call node.

        Args
        -----
        node : ast.Call
            The call AST node whose target name should be extracted.

        Returns
        --------
        str
            A best-effort string representing the function being called (e.g., "foo" or "obj.method").

        Notes
        -------
        - For attribute calls returns "<owner>.<attr>" if possible; otherwise returns a fallback "func".
        """
        if isinstance(node.func, ast.Name):
            return node.func.id
        elif isinstance(node.func, ast.Attribute):
            return getattr(node.func.value, 'id', 'obj') + "." + node.func.attr
        return "func"

def generate_function_flowchart(func_obj) -> str:
    """
    Description
    ----------
    Generate a Mermaid data flow diagram for a single Python function using AST analysis.

    Args
    -----
    func_obj : callable
        The function object to analyze.

    Returns
    --------
    str
        Mermaid markup representing the function's data flow graph (empty string on failure).

    Notes
    -------
    - Uses inspect.getsource and ast.parse to build an AST, then visits it with AdvancedFlowVisitor.
    - Returns an empty string if source extraction or parsing fails (e.g., builtins or dynamic functions).
    """
    try:
        source = inspect.getsource(func_obj)
        source = inspect.cleandoc(source)
        tree = ast.parse(source)
    except (OSError, TypeError, IndentationError, SyntaxError):
        return ""

    visitor = AdvancedFlowVisitor()
    visitor.visit(tree)

    if not visitor.nodes:
        return ""

        # Build Mermaid
    lines = ["flowchart TD"] # Use top-down layout, suitable for showing flow
    
    # Style definitions [Modified for Mental Model]
    # Input: Blue (Data Entry)
    lines.append("    classDef input fill:#e3f2fd,stroke:#1565c0,stroke-width:2px,rx:5,ry:5;")
    # Process: Gray/White (Standard Glue Code)
    lines.append("    classDef process fill:#fff,stroke:#bdbdbd,stroke-width:1px;")
    # Core Process: Yellow/Gold (The "Magic" happens here - AI/Math/IO)
    lines.append("    classDef core_process fill:#fff9c4,stroke:#fbc02d,stroke-width:2px,rx:5,ry:5;") 
    # Decision: Purple (Logic Branching)
    lines.append("    classDef decision fill:#f3e5f5,stroke:#7b1fa2,stroke-width:1px,rx:5,ry:5,stroke-dasharray: 5 5;")
    # Output: Green (Result)
    lines.append("    classDef output fill:#e8f5e9,stroke:#2e7d32,stroke-width:2px,rx:5,ry:5;")
    
    # Draw nodes
    for node in visitor.nodes:
        # 1. Sanitize ID
        safe_node_id = sanitize_id(node.id) 
        
        # 2. Escape special characters in Label
        # Replace newlines with <br/> for HTML labels, and escape quotes
        safe_label = node.label.replace('"', "'").replace('\n', '<br/>')
        
        shape_start, shape_end = "[", "]"
        if node.node_type == "input": shape_start, shape_end = "([", "])"
        if node.node_type == "output": shape_start, shape_end = "([", "])"
        if node.node_type == "decision": shape_start, shape_end = "{{", "}}" # Rhombus for decision
        
        # Use quotes around Label to ensure special characters (like spaces, =) are displayed correctly
        lines.append(f'    {safe_node_id}{shape_start}"{safe_label}"{shape_end}:::{node.node_type}')
        
        for source_id, var_name in node.edges_in:
            safe_source_id = sanitize_id(source_id)
            # Edge labels also need to be sanitized to remove characters that might break syntax
            safe_var = var_name.replace('"', "'").replace('|', '/').replace('\n', '')
            
            # Use dotted line for decision flows to distinguish them
            arrow = "-.->" if "Decision" in source_id else "-->"
            # Fix Mermaid syntax error
            # Incorrect: A --> "label" B
            # Correct: A -->|label| B
            if safe_var:
                lines.append(f'    {safe_source_id} {arrow}|"{safe_var}"| {safe_node_id}')
            else:
                lines.append(f'    {safe_source_id} {arrow} {safe_node_id}')

    return "\n".join(lines)

# --- Helper: 2. Global Call Graph Analysis (Macro) ---

class GlobalCallGraphVisitor(ast.NodeVisitor):
    """
    Description
    ----------
    Analyze the entire module's AST to build a call graph between functions.
    Traverses the AST to find function definitions and function calls, linking callers to callees.
    """
    def __init__(self, known_functions: Set[str]):
        """
        Description
        ----------
        Initialize the visitor with a set of known internal functions.

        Args
        -----
        known_functions : Set[str]
            A set of function names defined within the library being analyzed. Used to filter out built-in or external calls.

        Returns
        --------
        None

        Notes
        -------
        - Initializes `self.calls` to store the graph edges.
        - Initializes `self.dependency_map` for closure calculation.
        - Sets `self.current_function` to "Main_Script" to capture top-level calls.
        """
        self.known_functions = known_functions # Set of all function names defined in the library
        self.calls = [] # List of (caller, callee, arg_names)
        self.current_function = "Main_Script" # Default to top-level script
        # New: adjacency list for dependency closure calculation
        self.dependency_map = defaultdict(set) # caller -> set(callees)

    def visit_FunctionDef(self, node):
        """
        Description
        ----------
        Handle function definitions to track the current caller context.

        Args
        -----
        node : ast.FunctionDef
            The function definition AST node.

        Returns
        --------
        None

        Notes
        -------
        - Updates `self.current_function` to the name of the function being visited.
        - Restores the previous function name after visiting the body (handling nested functions).
        """
        prev_function = self.current_function
        self.current_function = node.name
        self.generic_visit(node)
        self.current_function = prev_function

    def visit_Call(self, node):    
        """
        Description
        ----------
        Handle function calls to record dependencies between the current function and the callee.

        Args
        -----
        node : ast.Call
            The function call AST node.

        Returns
        --------
        None

        Notes
        -------
        - Extracts the callee name (handling simple names and attributes like `self.method`).
        - Extracts argument names for visualization.
        - Filters calls: only records calls to functions in `known_functions` or `self.*` calls.
        - Updates `self.calls` and `self.dependency_map`.
        """
        # Extract the name of the called function
        callee_name = ""
        if isinstance(node.func, ast.Name):
            callee_name = node.func.id
        elif isinstance(node.func, ast.Attribute):
            # Handle self.method() or module.func()
            callee_name = node.func.attr
        
        if callee_name:
            # Extract argument names (for data flow visualization)
            args = []
            for arg in node.args:
                if isinstance(arg, ast.Name):
                    args.append(arg.id)
            
            # Only record calls to functions defined in our library (to avoid including built-ins like print, len)
            # Or if it's a self.xxx call, we also record it (assuming it's an internal class call)
            if callee_name in self.known_functions or (isinstance(node.func, ast.Attribute) and isinstance(node.func.value, ast.Name) and node.func.value.id == 'self'):
                self.calls.append((self.current_function, callee_name, args))
                self.dependency_map[self.current_function].add(callee_name)
        
        self.generic_visit(node)

# --- Helper: 3. Entry Analysis (Navigator: How to Drive) ---
class EntryAnalysisVisitor(ast.NodeVisitor):
    """
    Description
    ----------
    Analyze entry points in the code, especially the use of argparse.
    """
    def __init__(self):
        """
        Description
        ----------
        Initialize visitor state for entry point analysis.

        Args
        -----
        None

        Returns
        --------
        None

        Notes
        -------
        - Initializes `self.has_main_block` to False.
        - Initializes `self.args` to store discovered CLI arguments.
        """
        self.has_main_block = False
        self.args = [] # List of (arg_name, help_text)

    def visit_If(self, node):
        """
        Description
        ----------
        Detect `if __name__ == "__main__":` blocks to identify script entry points.

        Args
        -----
        node : ast.If
            The If statement AST node.

        Returns
        --------
        None

        Notes
        -------
        - Sets `self.has_main_block` to True if the standard main guard is found.
        """
        # Detect if __name__ == "__main__":
        try:
            if (isinstance(node.test, ast.Compare) and 
                isinstance(node.test.left, ast.Name) and 
                node.test.left.id == "__name__" and 
                isinstance(node.test.comparators[0], ast.Constant) and 
                node.test.comparators[0].value == "__main__"):
                self.has_main_block = True
        except: pass
        self.generic_visit(node)

    def visit_Call(self, node):
        """
        Description
        ----------
        Detect `parser.add_argument(...)` calls to extract CLI argument definitions.

        Args
        -----
        node : ast.Call
            The function call AST node.

        Returns
        --------
        None

        Notes
        -------
        - Extracts the argument name (e.g., '--input') and help text.
        - Appends found arguments to `self.args`.
        """
        # Detect parser.add_argument(...)
        if isinstance(node.func, ast.Attribute) and node.func.attr == 'add_argument':
            arg_name = "Unknown"
            help_text = ""
            
            # Extract argument name (usually the first argument)
            if node.args:
                if isinstance(node.args[0], ast.Constant):
                    arg_name = node.args[0].value
            
            # Extract help information
            for kw in node.keywords:
                if kw.arg == 'help' and isinstance(kw.value, ast.Constant):
                    help_text = kw.value.value
            
            self.args.append((arg_name, help_text))
        
        self.generic_visit(node)

# --- Helper: 4. Dependency Closure Calculation (Navigator: Extraction Guide) ---
def get_dependency_closure(target_func: str, dependency_map: Dict[str, Set[str]]) -> Set[str]:
    """
    Description
    ----------
    Calculate the dependency closure of the target function (i.e., all other functions required to run it).

    Args
    -----
    target_func : str
        The name of the function to analyze.
    dependency_map : Dict[str, Set[str]]
        Adjacency list representing the call graph (caller -> callees).

    Returns
    --------
    Set[str]
        A set of function names representing the transitive closure of dependencies.

    Notes
    -------
    - Uses Breadth-First Search (BFS) to traverse the dependency graph.
    """
    closure = set()
    queue = deque([target_func])
    visited = set()

    while queue:
        current = queue.popleft()
        if current in visited: continue
        visited.add(current)
        closure.add(current)

        if current in dependency_map:
            for dep in dependency_map[current]:
                if dep not in visited:
                    queue.append(dep)
    return closure

# --- Helper: 5. Module Classifier (Navigator: Architecture) ---
def classify_module(module_obj) -> str:
    """
    Description
    ----------
    Infer the role of a module based on the external libraries it imports or keywords in source.

    Args
    -----
    module_obj : module
        The module object to classify.

    Returns
    --------
    str
        A string category (e.g., "Model / AI", "Web / API", "Utility / Core").

    Notes
    -------
    - Scans source code for specific library imports (torch, flask, pandas, etc.).
    - Returns "Unknown" if source cannot be retrieved.
    """
    try:
        source = inspect.getsource(module_obj)
    except:
        return "Unknown"
    
    # Simple keyword matching
    if "torch" in source or "tensorflow" in source or "keras" in source:
        return "Model / AI"
    if "flask" in source or "django" in source or "fastapi" in source:
        return "Web / API"
    if "pandas" in source or "numpy" in source or "csv" in source:
        return "Data Processing"
    if "matplotlib" in source or "seaborn" in source or "plotly" in source:
        return "Visualization"
    if "argparse" in source or "click" in source:
        return "Interface / CLI"
    
    return "Utility / Core"

# --- Helper: 6. Safe Package Walker ---
def safe_walk_packages(path, prefix):
    """
    Description
    ----------
    Recursively walk packages like pkgutil.walk_packages, but aggressively skip
    test directories to avoid side-effects (like FileNotFoundError) during import.

    Args
    -----
    path : List[str]
        Paths to search for modules.
    prefix : str
        Prefix to output module names with.

    Yields
    -------
    Tuple[Any, str, bool]
        (importer, module_name, is_package)
    """
    skip_keywords = {'tests', 'test', 'conftest', 'examples', 'docs', 'setup', 'test_suite'}
    
    # pkgutil.iter_modules is safe, it scans directories without importing
    for importer, name, ispkg in pkgutil.iter_modules(path, prefix):
        # Check leaf name (e.g., 'metapredict.tests' -> 'tests')
        leaf_name = name.split('.')[-1]
        if leaf_name in skip_keywords:
            continue
            
        yield importer, name, ispkg
        
        if ispkg:
            try:
                # We must import the package to get its __path__ for recursion
                # Since we filtered 'tests' above, we hope this import is safe
                mod = importlib.import_module(name)
                if hasattr(mod, "__path__"):
                    yield from safe_walk_packages(mod.__path__, name + ".")
            except (KeyboardInterrupt, SystemExit):
                raise
            except:
                # Ignore import errors during recursion
                pass


def generate_global_call_graph(modules: List[Any], library_name: str) -> Tuple[str, Dict[str, Set[str]]]:
    """
    Description
    ----------
    Generate a global function call graph and return the dependency mapping.

    Args
    -----
    modules : List[Any]
        List of module objects to analyze.
    library_name : str
        The name of the library being inspected (used for filtering).

    Returns
    --------
    Tuple[str, Dict[str, Set[str]]]
        - str: Mermaid graph definition string.
        - Dict[str, Set[str]]: Dependency map (caller -> set of callees).

    Notes
    -------
    - 1. Collects all defined function names to create a whitelist.
    - 2. Traverses all source code using `GlobalCallGraphVisitor`.
    - 3. Constructs a Mermaid graph string with proper styling and ID sanitization.
    """
    # 1. Collect all defined function names (whitelist)
    known_functions = set()
    for mod in modules:
        for name, obj in inspect.getmembers(mod):
            if inspect.isfunction(obj) or inspect.ismethod(obj):
                known_functions.add(name)
            elif inspect.isclass(obj):
                for m_name, m_obj in inspect.getmembers(obj):
                    if inspect.isfunction(m_obj) or inspect.ismethod(m_obj):
                        known_functions.add(m_name)

    # 2. Traverse all source code for AST analysis
    visitor = GlobalCallGraphVisitor(known_functions)
    
    for mod in modules:
        try:
            source = inspect.getsource(mod)
            tree = ast.parse(source)
            visitor.visit(tree)
        except Exception:
            continue

    if not visitor.calls:
        return "", {}

    # 3. Construct Mermaid graph
    lines = ["graph TD"]
    lines.append("    classDef main fill:#f9f,stroke:#333,stroke-width:2px;")
    lines.append("    classDef func fill:#fff,stroke:#333,stroke-width:1px;")
    
    edges = set()
    
    # --- Create ID mapping to avoid keyword conflicts (e.g. function named 'graph') ---
    # Collect all unique function names involved in calls
    all_funcs = set()
    for caller, callee, _ in visitor.calls:
        all_funcs.add(caller)
        all_funcs.add(callee)
    
    # Map each function name to a safe ID like f_0, f_1...
    # This completely decouples the Mermaid ID from the function name
    id_map = {name: f"f_{i}" for i, name in enumerate(sorted(list(all_funcs)))}
    
    for caller, callee, args in visitor.calls:
        # Ignore recursive calls
        if caller == callee: continue
        
        # Use mapped IDs instead of sanitize_id(name)
        caller_id = id_map[caller]
        callee_id = id_map[callee]

        # Format edge
        edge_label = ""
        if args:
            # Truncate long argument lists to prevent graph explosion
            arg_str = '<br>'.join(args) 
            # if len(arg_str) > 20:
                # arg_str = arg_str[:17] + "..." 
            # Remove characters that may break Mermaid syntax
            arg_str = arg_str.replace('"', "'").replace('|', '/')
            edge_label = f"|{arg_str}|"
        
        # Use ID[Label] format
        # This way, the ID is safe (e.g., f_5), and the label can contain keywords like "graph"
        edge_str = f'    {caller_id}["{caller}"] -->{edge_label} {callee_id}["{callee}"]'
        
        if edge_str not in edges:
            edges.add(edge_str)
            lines.append(edge_str)
            
            # Apply styles to IDs
            if caller == "main" or caller == "Main_Script":
                lines.append(f"    class {caller_id} main;")
            else:
                lines.append(f"    class {caller_id} func;")
            lines.append(f"    class {callee_id} func;")

    return "\n".join(lines), visitor.dependency_map


def convert_md_to_html(md_content: str, title: str) -> str:
    """
    Description
    ----------
    Convert Markdown content to HTML with Mermaid rendering support.
    Improved Markdown parsing:
    - Uses Regex for inline formatting (bold, code) to fix partial line styling.
    - Implements a simple table parser to render HTML tables instead of raw text.
    - Fixes tag escaping issues for Mermaid.

    Args
    -----
    md_content : str
        The raw Markdown content string.
    title : str
        The title for the generated HTML page.

    Returns
    --------
    str
        A complete HTML string containing the rendered content and Mermaid.js integration.

    Notes
    -------
    - Splits content by code blocks to handle escaping differently for text vs code.
    - Manually restores specific HTML tags (details, summary) used in the report.
    - Escapes Mermaid diagram code to prevent browser parsing errors before Mermaid.js runs.
    """
    # Helper: Inline formatting using Regex
    def format_inline(text):
        # Bold: **text** -> <b>text</b>
        text = re.sub(r'\*\*(.*?)\*\*', r'<b>\1</b>', text)
        # Inline code: `text` -> <code>text</code>
        text = re.sub(r'`(.*?)`', r'<code>\1</code>', text)
        return text

    # Helper: Render a list of pipe-separated lines as an HTML table
    def render_table(buffer):
        if not buffer: return ""
        res = ["<table>"]
        # Header (Assume first line is header)
        headers = [h.strip() for h in buffer[0].strip().strip('|').split('|')]
        res.append("<thead><tr>" + "".join(f"<th>{format_inline(h)}</th>" for h in headers) + "</tr></thead>")
        res.append("<tbody>")
        # Body (Skip second line if it's a separator like |---|)
        start_idx = 2 if len(buffer) > 1 and '---' in buffer[1] else 1
        for row in buffer[start_idx:]:
            # Skip empty rows
            if not row.strip(): continue
            cells = [c.strip() for c in row.strip().strip('|').split('|')]
            res.append("<tr>" + "".join(f"<td>{format_inline(c)}</td>" for c in cells) + "</tr>")
        res.append("</tbody></table>")
        return "".join(res)

    parts = md_content.split("```")
    final_html_body = []
    
    for i, part in enumerate(parts):
        if i % 2 == 0:
            # === Plain text block ===
            text = html.escape(part)
            
            # Restore specific HTML tags
            text = text.replace("&lt;details&gt;", "<details>")
            text = text.replace("&lt;/details&gt;", "</details>")
            text = text.replace("&lt;summary&gt;", "<summary>")
            text = text.replace("&lt;/summary&gt;", "</summary>")
            
            lines = text.split('\n')
            formatted_lines = []
            table_buffer = [] # Buffer to hold table lines
            
            for line in lines:
                # --- Table Detection ---
                if line.strip().startswith('|'):
                    table_buffer.append(line)
                    continue
                
                # If we were in a table but now hit a non-table line, flush the table
                if table_buffer:
                    formatted_lines.append(render_table(table_buffer))
                    table_buffer = []
                
                # --- Normal Line Processing ---
                # Apply inline formatting (Bold, Code)
                line_fmt = format_inline(line)
                
                if line.startswith('# '): formatted_lines.append(f"<h1>{line_fmt[2:]}</h1>")
                elif line.startswith('## '): formatted_lines.append(f"<h2>{line_fmt[3:]}</h2>")
                elif line.startswith('### '): formatted_lines.append(f"<h3>{line_fmt[4:]}</h3>")
                elif line.startswith('#### '): formatted_lines.append(f"<h4>{line_fmt[5:]}</h4>")
                # Note: html.escape converts > to &gt;
                elif line.startswith('&gt; '): formatted_lines.append(f"<blockquote>{line_fmt[5:]}</blockquote>")
                elif line.strip() == "": formatted_lines.append("<br>")
                else: formatted_lines.append(f"{line_fmt}<br>")
            
            # Flush any remaining table at the end of the block
            if table_buffer:
                formatted_lines.append(render_table(table_buffer))
            
            final_html_body.append("\n".join(formatted_lines))
        else:
            # === Code block ===
            if part.startswith("mermaid"):
                # Mermaid diagram
                graph_code = part[7:].strip()
                escaped_code = html.escape(graph_code)
                final_html_body.append(f'<div class="mermaid" style="overflow-x: auto;">\n{escaped_code}\n</div>')
            else:
                # Plain code
                lang = part.split('\n')[0]
                code = part[len(lang):].strip()
                escaped_code = html.escape(code)
                final_html_body.append(f'<pre style="background:#f4f4f4; padding:10px; border-radius:5px;"><code>{escaped_code}</code></pre>')

    body_str = "\n".join(final_html_body)

    html_template = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif; line-height: 1.6; padding: 20px; max-width: 1200px; margin: 0 auto; color: #333; }}
        h1, h2, h3 {{ color: #24292e; border-bottom: 1px solid #eaecef; padding-bottom: .3em; }}
        code {{ background-color: #f6f8fa; padding: 0.2em 0.4em; border-radius: 3px; font-family: monospace; }}
        pre {{ background-color: #f6f8fa; padding: 16px; overflow: auto; border-radius: 6px; }}
        blockquote {{ border-left: 4px solid #dfe2e5; color: #6a737d; padding-left: 1em; margin-left: 0; }}
        
        /* Table Styles */
        table {{ border-collapse: collapse; width: 100%; margin-bottom: 16px; display: block; overflow-x: auto; }}
        th, td {{ border: 1px solid #dfe2e5; padding: 6px 13px; }}
        th {{ background-color: #f6f8fa; font-weight: 600; }}
        tr:nth-child(2n) {{ background-color: #f6f8fa; }}
        
        .mermaid {{ margin: 20px 0; text-align: center; }}
        details {{ margin-bottom: 10px; border: 1px solid #e1e4e8; border-radius: 6px; padding: 8px; }}
        summary {{ cursor: pointer; font-weight: bold; outline: none; }}
        
        .edgeLabel {{
            font-size: 11px !important;
            background-color: rgba(255, 255, 255, 0.9) !important;
            padding: 2px !important;
            border-radius: 4px;
        }}
    </style>
</head>
<body>
    {body_str}

    <!-- Import Mermaid.js -->
    <script type="module">
        import mermaid from 'https://cdn.jsdelivr.net/npm/mermaid@10/dist/mermaid.esm.min.mjs';
        mermaid.initialize({{ 
            startOnLoad: true, 
            maxTextSize: 1000000,  
            maxEdges: 100000, 
            theme: 'default',
            flowchart: {{ 
                useMaxWidth: false, 
                htmlLabels: true,
                rankSpacing: 150, 
                nodeSpacing: 100,
                curve: 'basis' 
            }} 
        }});
    </script>
</body>
</html>
    """
    return html_template


def inspect_library(
    library_name: str,
    output_path: Optional[str] = None,
    include_private: bool = False,
    include_imported: bool = False,
    limit_api_recommendations: Optional[int] = 20,
    limit_code_snippets: Optional[int] = 20,
    limit_snippet_args: Optional[int] = 3,
    limit_external_libs: Optional[int] = 20,
    limit_pagerank: Optional[int] = 20,
    limit_dependency_graph: Optional[int] = 100,
    limit_inheritance_graph: Optional[int] = 100,
    limit_extraction_guide: Optional[int] = 20,
    run_all_candidates: bool = True,
    parent_resolution: Optional[object] = None,
):
    """
    Description
    ----------
    Main entry point to inspect a Python library.
    Performs dynamic import, AST analysis, dependency graph construction, and report generation.
    
    [Updated] Now supports CLI commands, local files, and PyPI package names via TargetResolver.

    Args
    -----
    library_name : str
        The importable name of the library to inspect (e.g., "pandas", "requests").
    output_path : Optional[str]
        Path to save the generated Markdown report. If None, prints to stdout.
    include_private : bool
        If True, includes private members (starting with _) in the report.
    include_imported : bool
        If True, includes members imported from other libraries in the API list.
    limit_api_recommendations : int
        Limit for the number of recommended API entry points.
    limit_code_snippets : int
        Limit for the number of auto-generated code snippets.
    limit_snippet_args : int
        Threshold for argument count to switch to multi-line formatting in snippets.
    limit_external_libs : int
        Limit for the number of top external dependencies shown.
    limit_pagerank : int
        Limit for the number of top modules by PageRank shown.
    limit_dependency_graph : int
        Limit for the number of nodes in the dependency graph.
    limit_inheritance_graph : int
        Limit for the number of classes in the inheritance graph.
    limit_extraction_guide : int
        Limit for the number of functions in the extraction guide.
    run_all_candidates : bool
        If True, attempts to run inspection on all resolved import candidates (for CLI commands or files).

    Returns
    --------
    None

    Notes
    -------
    - Phase 0: Navigator - Analyzes entry points (CLI or API) and suggests usage.
    - Phase 1: Network Analysis - Builds dependency graphs (internal & external).
    - Phase 2: API Inspection - Lists functions/classes and generates logic flowcharts.
    - Generates both Markdown and HTML reports if output_path is provided.
    """
    # --- 0. Resolve Target ---
    resolver = TargetResolver()
    resolved = resolver.resolve(library_name)
    
    real_lib_name = resolved.import_name
    print(f"🔍 Target Resolution: Input='{library_name}' -> Import='{real_lib_name}' (Type: {resolved.target_type})")

    # --- 1. Dynamically import the main library (with sys.argv protection) ---
    _old_argv = sys.argv
    sys.argv = [sys.argv[0]]

    submodules = []
    main_module = None

    # Build ordered import candidates to maximize chance of successful import
    import_candidates = []  # list of (candidate_name, hint)
    file_candidate = None

    # Prefer the resolved import name
    if resolved.import_name:
        import_candidates.append((resolved.import_name, 'resolved'))

    # If this was a CLI command resolution, include the concrete entry module (ep.module)
    # but append it after the resolved package/module candidates so that
    # single-run mode prefers the package name first.
    if resolved.target_type == 'cli_command' and getattr(resolved, 'cli_info', None):
        try:
            ep_module = resolved.cli_info.get('module')
            if ep_module and ep_module not in [c for c, _ in import_candidates]:
                import_candidates.append((ep_module, 'entry_point_module'))
        except Exception:
            pass

    # If target is a file, remember file path and attempt to construct a dotted module name
    if resolved.target_type == 'file' and resolved.file_path:
        file_candidate = resolved.file_path
        abs_fp = os.path.abspath(resolved.file_path)
        # Candidate: module name relative to repo root if possible
        try:
            # repo root guess: look for pyproject.toml or setup.py upwards from file
            def find_package_root(start_path):
                cur = os.path.abspath(start_path)
                last = None
                while True:
                    if os.path.exists(os.path.join(cur, 'pyproject.toml')) or os.path.exists(os.path.join(cur, 'setup.py')):
                        return cur
                    if os.path.exists(os.path.join(cur, '__init__.py')):
                        # package dir found; use its parent as root
                        return os.path.dirname(cur)
                    parent = os.path.dirname(cur)
                    if parent == cur or parent == last:
                        return None
                    last = cur
                    cur = parent

            pkg_root = find_package_root(os.path.dirname(abs_fp))
            if pkg_root:
                mod_rel = os.path.relpath(abs_fp, pkg_root)
                mod_name = os.path.splitext(mod_rel)[0].replace(os.sep, '.')
                if mod_name not in [c for c, _ in import_candidates]:
                    import_candidates.insert(0, (mod_name, 'constructed_from_file'))
                # also ensure pkg_root is considered later when importing
            else:
                # fallback: use file stem as module
                mod_name = os.path.splitext(os.path.basename(abs_fp))[0]
                if mod_name not in [c for c, _ in import_candidates]:
                    import_candidates.append((mod_name, 'file_stem'))
        except Exception:
            pass

    # Ensure uniqueness while preserving order
    seen = set()
    uniq_imports = []
    for name, hint in import_candidates:
        if name in seen: continue
        seen.add(name)
        uniq_imports.append((name, hint))

    # If requested, run inspect on each candidate separately to maximize coverage
    if run_all_candidates:
        candidates = [name for name, _ in uniq_imports]
        if file_candidate:
            candidates.append(file_candidate)

        if output_path:
            base, ext = os.path.splitext(output_path)
            if not ext:
                ext = '.md'
        else:
            base = None

        for i, cand in enumerate(candidates, start=1):
            if base:
                safe = re.sub(r'[^a-zA-Z0-9_.-]', '_', os.path.basename(str(cand)))
                out = f"{base}_candidate_{i}_{safe}{ext}"
            else:
                out = None
            print(f"\n🔁 Running candidate {i}/{len(candidates)}: {cand} -> output={out}")
            try:
                inspect_library(
                    str(cand),
                    output_path=out,
                    include_private=include_private,
                    include_imported=include_imported,
                    limit_api_recommendations=limit_api_recommendations,
                    limit_code_snippets=limit_code_snippets,
                    limit_snippet_args=limit_snippet_args,
                    limit_external_libs=limit_external_libs,
                    limit_pagerank=limit_pagerank,
                    limit_dependency_graph=limit_dependency_graph,
                    limit_inheritance_graph=limit_inheritance_graph,
                    limit_extraction_guide=limit_extraction_guide,
                    run_all_candidates=False,
                    parent_resolution=resolved,
                )
            except Exception as e:
                print(f"⚠️ Candidate {cand} failed: {e}")
        return

    try:
        # Try each import candidate until one succeeds
        imported_ok = False
        tried_candidates = []  # records of (name, hint, status, error)
        used_candidate_name = None
        used_candidate_hint = None
        for cand_name, hint in uniq_imports:
            try:
                # If file_candidate exists and hint suggests constructed_from_file, add its root to sys.path
                if file_candidate and hint in ('constructed_from_file', 'file_stem'):
                    # attempt to find package root again and insert into sys.path
                    cur_dir = os.path.dirname(os.path.abspath(file_candidate))
                    # walk up up to 4 levels and add to sys.path for best chance
                    temp = cur_dir
                    for _ in range(4):
                        if temp not in sys.path:
                            sys.path.insert(0, temp)
                        temp = os.path.dirname(temp)

                main_module = importlib.import_module(cand_name)
                submodules.append(main_module)
                imported_ok = True
                real_lib_name = cand_name
                used_candidate_name = cand_name
                used_candidate_hint = hint
                tried_candidates.append((cand_name, hint, 'success', ''))
                break
            except ImportError as e:
                tried_candidates.append((cand_name, hint, 'failed', str(e)))
                # try next candidate
                continue
            except Exception as e:
                tried_candidates.append((cand_name, hint, 'failed', str(e)))
                continue

        if not imported_ok:
            # Last resort: if resolved is file and file path exists, try adding its directory to sys.path and import by stem
            if file_candidate and os.path.isfile(file_candidate):
                fp_dir = os.path.dirname(os.path.abspath(file_candidate))
                if fp_dir not in sys.path:
                    sys.path.insert(0, fp_dir)
                try:
                    stem = os.path.splitext(os.path.basename(file_candidate))[0]
                    main_module = importlib.import_module(stem)
                    submodules.append(main_module)
                    imported_ok = True
                    real_lib_name = stem
                    used_candidate_name = stem
                    used_candidate_hint = 'file_stem_attempt'
                    tried_candidates.append((stem, 'file_stem_attempt', 'success', ''))
                except Exception as e:
                    # Still failed; record and abort import-phase (we may still attempt file-level analysis downstream)
                    tried_candidates.append((file_candidate, 'file_path_attempt', 'failed', str(e)))
                    print(f"❌ Error: Could not import library '{resolved.import_name if resolved.import_name else library_name}'. Reason: {e}")
                    if resolved.target_type == 'cli_command':
                        print(f"   (Tip: The CLI command '{library_name}' maps to module '{resolved.import_name}', which seems missing or broken.)")
                    # fall through to allow file-level analysis if implemented

        if not imported_ok:
            # If no import succeeded, attempt a file-level AST fallback (if we have a file path)
            if file_candidate and os.path.isfile(file_candidate):
                try:
                    src = open(file_candidate, 'r', encoding='utf-8').read()
                    tree = ast.parse(src)

                    fb_lines = []
                    fb_lines.append(f"# File-level Documentation for `{os.path.basename(file_candidate)}`")
                    fb_lines.append(f"> **Note:** Import failed; performing AST-only file analysis for `{file_candidate}`.")
                    fb_lines.append(f"**File Path:** `{file_candidate}`\n")
                    run_mode = 'AST-only fallback'

                    # Module docstring
                    file_doc = ast.get_docstring(tree)
                    if file_doc:
                        fb_lines.append("## Module Docstring")
                        fb_lines.append(f"```text\n{file_doc}\n```\n")

                    # Collect imports
                    import_counter = Counter()
                    for node in ast.walk(tree):
                        if isinstance(node, ast.Import):
                            for n in node.names:
                                import_counter[n.name] += 1
                        elif isinstance(node, ast.ImportFrom):
                            module = node.module if node.module else "<relative>"
                            import_counter[module] += 1

                    if import_counter:
                        fb_lines.append("## Top Imports Detected")
                        fb_lines.append("| Module | Count |")
                        fb_lines.append("| :--- | :---: |")
                        for mod, cnt in import_counter.most_common():
                            fb_lines.append(f"| `{mod}` | {cnt} |")
                        fb_lines.append("")

                    # Top-level functions and classes
                    funcs = []
                    classes = []
                    for node in tree.body:
                        if isinstance(node, ast.FunctionDef):
                            funcs.append(node)
                        elif isinstance(node, ast.ClassDef):
                            classes.append(node)

                    if funcs:
                        fb_lines.append("## Top-level Functions")
                        for fn in funcs:
                            # signature approximation
                            args = []
                            for a in fn.args.args:
                                args.append(a.arg)
                            sig = f"({', '.join(args)})"
                            fb_lines.append(f"- `{fn.name}{sig}`")
                        fb_lines.append("")

                        # Include small snippets for each function
                        for fn in funcs:
                            fb_lines.append(f"### `{fn.name}`")
                            try:
                                src_snip = ast.get_source_segment(src, fn)
                                if not src_snip:
                                    src_snip = "# source unavailable"
                            except Exception:
                                src_snip = "# source unavailable"
                            fb_lines.append("```python")
                            fb_lines.append(src_snip)
                            fb_lines.append("```")
                            fb_lines.append("")

                    if classes:
                        fb_lines.append("## Top-level Classes")
                        for cl in classes:
                            fb_lines.append(f"- `{cl.name}`")
                        fb_lines.append("")

                    # Entry detection
                    entry_vis = EntryAnalysisVisitor()
                    try:
                        entry_vis.visit(tree)
                        if entry_vis.has_main_block:
                            fb_lines.append("- ✅ **Script Entry Point**: Contains `if __name__ == '__main__'` block.")
                    except Exception:
                        pass

                    md_content = "\n".join(fb_lines)
                    # Write markdown output
                    if output_path:
                        out = output_path
                        # ensure directory exists
                        od = os.path.dirname(out)
                        if od and not os.path.exists(od):
                            os.makedirs(od, exist_ok=True)
                        with open(out, 'w', encoding='utf-8') as f:
                            f.write(md_content)
                        print(f"✅ Markdown report saved to: {out}")
                        # Also write HTML version
                        try:
                            html_out = os.path.splitext(out)[0] + '.html'
                            html_content = convert_md_to_html(md_content, title=os.path.basename(out))
                            with open(html_out, 'w', encoding='utf-8') as hf:
                                hf.write(html_content)
                            print(f"📊 Interactive HTML report saved to: {html_out}")
                        except Exception:
                            pass
                    else:
                        print(md_content)

                    return
                except Exception as e:
                    print(f"❌ File-level AST analysis failed for '{file_candidate}': {e}")
                    return
            # No file fallback available; abort
            print(f"❌ Error: Unable to import any candidate modules for target '{library_name}'.")
            return

        print(f"🔍 Analyzing dependencies for '{real_lib_name}' (Network Analysis Phase)...")

        if hasattr(main_module, "__path__"):
            # Use safe_walk_packages instead of pkgutil.walk_packages
            # This prevents importing 'tests' packages that might have side effects (like missing files)
            for importer, modname, ispkg in safe_walk_packages(main_module.__path__, main_module.__name__ + "."):
                # Double check skip logic (though safe_walk_packages handles most)
                mod_parts = modname.split('.')
                if 'tests' in mod_parts or 'test' in mod_parts or 'conftest' in mod_parts:
                    continue

                try:
                    sub_mod = importlib.import_module(modname)
                    submodules.append(sub_mod)
                except (KeyboardInterrupt, SystemExit):
                    # Allow users to interrupt with Ctrl+C
                    raise
                except BaseException:
                    # Exceptions raised by pytest.skip() inherit from BaseException, not Exception
                    # This allows ignoring import errors caused by missing optional dependencies
                    continue
    finally:
        sys.argv = _old_argv

    lines = []
    lines.append(f"# Documentation for `{real_lib_name}`")

    # Discovery vs current resolution context
    if parent_resolution is not None:
        try:
            pr_type = getattr(parent_resolution, 'target_type', str(parent_resolution))
            pr_input = getattr(parent_resolution, 'original_input', None) or getattr(parent_resolution, 'import_name', None)
            lines.append(f"> **Discovered via:** `{pr_type}` (original input: `{pr_input}`)\n")
        except Exception:
            pass

    # Current resolution note
    try:
        cur_type = getattr(resolved, 'target_type', 'unknown')
    except Exception:
        cur_type = 'unknown'

    if cur_type == 'cli_command':
        lines.append(f"> **Note:** Analyzed via CLI command `{library_name}` (resolved as CLI entry).\n")
    elif cur_type == 'file':
        lines.append(f"> **Note:** Analyzed local file/package at `{resolved.original_input}`.\n")
    else:
        lines.append(f"> **Note:** Resolved as `{cur_type}` for `{real_lib_name}`.\n")

    lines.append(f"**File Path:** `{getattr(main_module, '__file__', 'Built-in/Unknown')}`\n")

    doc = None
    try:
        doc = inspect.getdoc(main_module)
    except Exception:
        doc = None
    if doc:
        lines.append("## Module Docstring")
        lines.append(f"```text\n{doc}\n```\n")

    # ===== Metadata & Diagnostics =====
    lines.append("## 🧾 Metadata & Diagnostics")

    # Used candidate
    if 'used_candidate_name' in locals() and used_candidate_name:
        lines.append(f"- **Used candidate:** `{used_candidate_name}` ({used_candidate_hint})")
    else:
        lines.append(f"- **Used candidate:** `None` (will run file-level AST fallback if available)")

    # Tried candidates table
    if 'tried_candidates' in locals() and tried_candidates:
        lines.append("### Tried Candidates")
        lines.append("| Candidate | Hint | Status | Error |")
        lines.append("| :--- | :--- | :--- | :--- |")
        for name, hint, status, err in tried_candidates:
            safe_err = (err.replace('\n', ' ')[:200] if err else '')
            lines.append(f"| `{name}` | {hint} | {status} | `{safe_err}` |")
        lines.append("")

    # Package / module version (best-effort)
    ver = None
    try:
        if imported_ok and main_module is not None:
            ver = getattr(main_module, '__version__', None)
        if not ver:
            try:
                ver = importlib_metadata.version(real_lib_name)
            except Exception:
                pass
    except Exception:
        ver = None
    if ver:
        lines.append(f"- **Package Version:** `{ver}`")

    # Run mode note
    if 'run_mode' in locals() and run_mode:
        lines.append(f"- **Run Mode:** {run_mode}")
    else:
        if imported_ok:
            lines.append("- **Run Mode:** dynamic import (module)")
        else:
            lines.append("- **Run Mode:** AST-only fallback (no import succeeded)")

    # Sys.path snapshot (helpful for diagnosing import issues)
    try:
        sp = sys.path[:10]  # show only top 10 entries
        lines.append("### sys.path (top 10 entries)")
        lines.append("```text")
        for p in sp:
            lines.append(p)
        lines.append("```\n")
    except Exception:
        pass

    # ==========================================
    # Phase 0: Navigator - How to Drive (Enhanced)
    # ==========================================
    lines.append("## 🚦 Navigator: How to Drive")
    lines.append("This section helps you understand how to run this library from the command line or entry points.")
    
    # --- System-level CLI Command Detection ---
    # Scan for console_scripts registered by this package
    detected_cli_commands = []
    try:
        # Try to find the distribution name (usually same as lib name, or mapped)
        dist_name = resolved.dist_name or real_lib_name
        
        # Get all console scripts
        eps = importlib_metadata.entry_points(group='console_scripts')
        
        # Filter for scripts belonging to this distribution OR this module
        for ep in eps:
            # Check if entry point points to the current module or submodules
            if ep.module.startswith(real_lib_name):
                detected_cli_commands.append(ep)
            # Or if we know the distribution name, check that
            elif hasattr(ep, 'dist') and ep.dist and ep.dist.name == dist_name:
                detected_cli_commands.append(ep)
    except Exception: pass

    if detected_cli_commands:
        lines.append("\n### 💻 Installed CLI Commands")
        lines.append("This library installs the following system commands (accessible from terminal):")
        lines.append("| Command | Entry Point (Function) |")
        lines.append("| :--- | :--- |")
        for ep in list(set(detected_cli_commands)): # Deduplicate
            lines.append(f"| `{ep.name}` | `{ep.value}` |")
        lines.append("")
        
        if resolved.target_type == 'cli_command':
             lines.append(f"- ✅ **Target Match**: You are analyzing the package backing the command `{resolved.original_input}`.")

    # ---  Script Entry Point Detection ---
    entry_visitor = EntryAnalysisVisitor()
    try:
        source = inspect.getsource(main_module)
        entry_visitor.visit(ast.parse(source))
    except: pass

    if entry_visitor.has_main_block:
        lines.append("- ✅ **Script Entry Point**: This module contains an `if __name__ == '__main__':` block, meaning it can be run directly.")
    elif not detected_cli_commands:
        lines.append("- ℹ️ **No Direct Entry Point**: This module seems to be a library intended for import, not direct execution.")
        
    # ---  Python API Usage (Inferred) ---
    lines.append("\n### 🐍 Python API Usage (Inferred)")
    lines.append("Since no CLI entry point was found, here are the likely **Python API entry points** for your script:")
    
    library_name = real_lib_name # Override for subsequent logic

    # 1. Pre-build call graph to assist judgment (partially run Phase 1 logic in advance)
    # We need to know which functions are "Top Level" (not called internally)
    temp_known_funcs = set()
    temp_calls = []
    if hasattr(main_module, "__path__"):
        # Simple scan of the main module and direct submodules
        scan_modules = [main_module] + submodules # Limit number to avoid being too slow submodules[:10]
    else:
        scan_modules = [main_module]

    for m in scan_modules:
        for n, o in inspect.getmembers(m):
            if inspect.isfunction(o) or inspect.ismethod(o): temp_known_funcs.add(n)
        
    temp_visitor = GlobalCallGraphVisitor(temp_known_funcs)
    for m in scan_modules:
        try: temp_visitor.visit(ast.parse(inspect.getsource(m)))
        except: pass
            
    # Calculate in-degree (number of times called)
    in_degree = Counter()
    for caller, callee, _ in temp_visitor.calls:
        in_degree[callee] += 1

    # 2. Find potential API entry points
    api_candidates = []
    candidates_source = []
        
    if hasattr(main_module, "__all__"):
        candidates_source = main_module.__all__
    else:
        candidates_source = [x for x in dir(main_module) if not x.startswith("_")]
            
    for name in candidates_source:
        try:
            obj = getattr(main_module, name)
            # Only focus on functions and classes
            if inspect.isfunction(obj):
                api_candidates.append((name, "Function", obj))
            elif inspect.isclass(obj):
                # Ignore exception classes
                if issubclass(obj, Exception): continue
                api_candidates.append((name, "Class", obj))
        except: continue
            
    # 3. Enhanced scoring mechanism
    def score_candidate(item):
        """
        Description
        ----------
        Calculate a relevance score for an API candidate to determine if it's a likely entry point.

        Args
        -----
        item : Tuple[str, str, Any]
            A tuple containing (name, kind, object), where kind is 'Function' or 'Class'.

        Returns
        --------
        int
            A score indicating the likelihood of being a main entry point (higher is better).

        Notes
        -------
        - Scoring factors:
            - Keywords in name (verbs like predict, run, load).
            - Class names (Model, Engine).
            - Topological usage (low in-degree suggests public API).
            - Code complexity (lines of code).
        """
        name, kind, obj = item
        score = 0
        name_lower = name.lower()
            
        # A. Keyword scoring (business logic priority)
        # Core verbs
        if any(v in name_lower for v in ['predict', 'calculate', 'analyze', 'solve', 'run', 'process', 'train', 'evaluate', 'generate']): score += 15
        # Entry verbs
        elif any(v in name_lower for v in ['main', 'start', 'init', 'load', 'create', 'make']): score += 10
        # Conversion/tool verbs
        elif any(v in name_lower for v in ['convert', 'parse', 'read', 'write', 'save', 'plot', 'show']): score += 5
            
        # B. Type scoring
        if kind == 'Class':
            # Model/core class bonus
            if any(n in name_lower for n in ['model', 'engine', 'client', 'api', 'runner', 'predictor']): score += 8
            
        # C. Topological scoring (key improvement)
        # If a function is public and rarely called internally (low in-degree), it is likely intended for external use
        if kind == 'Function':
            degree = in_degree.get(name, 0)
            if degree == 0: score += 10  # Pure top-level interface
            elif degree < 3: score += 5  # Low coupling interface
            else: score -= 5             # Heavily called internally, likely a low-level utility function
            
        # D. Complexity scoring (lines of code)
        try:
            lines_of_code = len(inspect.getsource(obj).splitlines())
            if lines_of_code > 50: score += 5 # Complex logic often indicates main interface
            elif lines_of_code < 3: score -= 5 # Wrappers usually have only one or two lines
        except: pass

        return score

    api_candidates.sort(key=score_candidate, reverse=True)
        
    # 4. Generate code snippets

    '''
    if api_candidates: 
        lines.append("```python")
        lines.append(f"import {library_name}")
        lines.append("")
        lines.append("# Likely entry points (Ranked by relevance):")
        # Show top 8 for better coverage
        for name, kind, obj in api_candidates: 
            try:
                sig = str(inspect.signature(obj))
            except: sig = "(...)"
                
            # Add a brief comment line (take the first line of the docstring)
            doc = inspect.getdoc(obj)
            doc_summary = f"  # {doc.splitlines()[0]}" if doc else ""
            # Truncate overly long comments
            if len(doc_summary) > 60: doc_summary = doc_summary[:57] + "..."
                
            lines.append(f"# {kind}: {library_name}.{name}{sig}{doc_summary}")
        lines.append("```")
    else:
        lines.append("_No obvious public API members detected._")
    '''

    # ---  User-friendly formatted output ---
    # ⚠️[ARG1]
    if api_candidates:
        lines.append(f"\n#### 🚀 Top {limit_api_recommendations} Recommended Entry Points")
        lines.append("| Type | API | Description |")
        lines.append("| :--- | :--- | :--- |")
            
        count = 0
        for name, kind, obj in api_candidates:
            # Filter out candidates with low scores (likely low-level utilities)
            if score_candidate((name, kind, obj)) < 0 and count > 5: continue
            if count >= limit_api_recommendations: break # Show up to limit_api_recommendations core APIs, default is 20, means that we only show top 20 APIs
                
            # 1. Extract signature and beautify parameters
            try:
                sig = inspect.signature(obj)
                params = []
                for p_name, p in sig.parameters.items():
                    if p.default == inspect.Parameter.empty:
                        # Required parameters in bold
                        params.append(f"**{p_name}**")
                    else:
                        # Optional parameters in regular font
                        params.append(f"{p_name}")
                    
                # Reassemble signature to avoid excessive length
                sig_str = f"({', '.join(params)})"
                # if len(sig_str) > 60: # If parameters are too long, truncate
                #    sig_str = f"({', '.join(params[:3])}, ...)"
            except: 
                sig_str = "(...)"

            # 2. Extract and clean documentation
            doc = inspect.getdoc(obj)
            doc_summary = doc.splitlines()[0] if doc else "No description."
            # if len(doc_summary) > 80: doc_summary = doc_summary[:77] + "..."
                
            # 3. Icon differentiation
            icon = "ƒ" if kind == "Function" else "C"
                
            lines.append(f"| `{icon}` | **{library_name}.{name}**{sig_str} | {doc_summary} |")
            count += 1
            
        lines.append("\n> **Note:** Bold parameters are required. Others are optional.")
            
        # --- New: Multi-dimensional code snippet generation ---
        lines.append("\n#### 🧩 Code Snippets (Auto-Generated)")
        lines.append("```python")
        lines.append(f"import {library_name}")
        lines.append("")
            
        # 1. Extract Top Functions (up to 20)
        # ⚠️[ARG2] 
        top_funcs = [x for x in api_candidates if x[1] == 'Function'][:limit_code_snippets]
        if top_funcs:
            lines.append(f"# --- Top {limit_code_snippets} Ranked Functions ---") # or we can just use {len(top_funcs)}
            for i, (name, _, obj) in enumerate(top_funcs):
                try:
                    sig = inspect.signature(obj)
                    args = []
                    for p_name, p in sig.parameters.items():
                        if p.default == inspect.Parameter.empty and p_name != 'self':
                            args.append(f"{p_name}=...") 
                        elif p.default != inspect.Parameter.empty:
                            # Optional parameters are commented out to indicate their presence
                            # args.append(f"{p_name}={p.default}") 
                            pass
                        
                    # If there are too many parameters, display them on multiple lines
                    # ⚠️[ARG3]
                    if len(args) > limit_snippet_args: # default is 3, switch to multi-line 
                        args_str = ",\n    ".join(args)
                        call_str = f"{library_name}.{name}(\n    {args_str}\n)"
                    else:
                        call_str = f"{library_name}.{name}({', '.join(args)})"
                            
                    lines.append(f"# {i+1}. {name}")
                    lines.append(f"result_{i+1} = {call_str}")
                    lines.append("")
                except: pass
        
        # 2. Extract Top Classes (up to 20)
        # ⚠️[ARG2]
        top_classes = [x for x in api_candidates if x[1] == 'Class'][:limit_code_snippets]
        if top_classes:
            lines.append(f"# --- Top {limit_code_snippets} Core Classes Initialization ---") # or we can just use {len(top_classes)}
            for i, (name, _, obj) in enumerate(top_classes):
                try:
                    # Attempt to get the signature of __init__
                    sig = inspect.signature(obj)
                    args = []
                    for p_name, p in sig.parameters.items():
                        if p.default == inspect.Parameter.empty and p_name != 'self':
                            args.append(f"{p_name}=...")
                        
                    call_str = f"{library_name}.{name}({', '.join(args)})"
                    instance_name = name.lower().replace("model", "_model").replace("class", "_obj")
                        
                    lines.append(f"# {i+1}. {name}")
                    lines.append(f"{instance_name} = {call_str}")
                    lines.append("")
                except: pass
            
        lines.append("```")

    else:
        lines.append("_No obvious public API members detected._")

    if entry_visitor.args:
        lines.append("\n### ⌨️ CLI Arguments (Detected)")
        lines.append("| Argument | Help Text |")
        lines.append("| :--- | :--- |")
        for arg, help_text in entry_visitor.args:
            lines.append(f"| `{arg}` | {help_text} |")
    else:
        lines.append("\n_No explicit `argparse` configuration detected in the main module._")
    lines.append("\n")


    # ==========================================
    # Phase 1: Network Construction & Analysis
    # ==========================================

    G = nx.DiGraph() if HAS_NETWORKX else None
    internal_modules_rank = Counter() 
    external_libs_rank = Counter()    
    dependency_graph = defaultdict(set) 
    module_roles = {} # module_name -> role

    for mod in submodules:
        current_mod_name = mod.__name__
        
        # Module role classification
        role = classify_module(mod)
        module_roles[current_mod_name] = role

        if HAS_NETWORKX:
            G.add_node(current_mod_name, type='internal', role=role)
        
        for name, obj in inspect.getmembers(mod):
            obj_module = getattr(obj, "__module__", None)
            if not obj_module: continue
            if obj_module == current_mod_name: continue
            dependency_graph[current_mod_name].add(obj_module)
            if obj_module.startswith(library_name):
                internal_modules_rank[obj_module] += 1
                if HAS_NETWORKX: G.add_edge(current_mod_name, obj_module)
            else:
                top_level_pkg = obj_module.split('.')[0]
                if top_level_pkg not in ['builtins', 'sys', 'os', 'typing']:
                    external_libs_rank[top_level_pkg] += 1
                    if HAS_NETWORKX:
                        G.add_node(top_level_pkg, type='external')
                        G.add_edge(current_mod_name, top_level_pkg)

    lines.append("## 📊 Network & Architecture Analysis")
    if not HAS_NETWORKX: lines.append("> ⚠️ `networkx` is not installed. Advanced metrics are disabled.\n")
    # ⚠️[ARG4]
    lines.append(f"### 🌍 Top {limit_external_libs} External Dependencies") # default is 20, now change to {limit_external_libs}
    if external_libs_rank:
        lines.append("| Library | Usage Count |")
        lines.append("| :--- | :--- |")
        for lib, count in external_libs_rank.most_common(limit_external_libs):
            lines.append(f"| **{lib}** | {count} |")
    else:
        lines.append("_No significant external dependencies._")
    lines.append("\n")
    
    sorted_pr = []
    if HAS_NETWORKX and len(G.nodes) > 0:
        lines.append("### 🕸️ Network Metrics (Advanced)")
        try:
            pagerank = nx.pagerank(G, alpha=0.85)
            sorted_pr = sorted(pagerank.items(), key=lambda x: x[1], reverse=True)
            # ⚠️[ARG5]
            lines.append(f"#### 👑 Top {limit_pagerank} Modules by PageRank (Authority)")
            lines.append("| Rank | Module | Score | Type | Role |")
            lines.append("| :--- | :--- | :--- | :--- | :--- |")
            for i, (node, score) in enumerate(sorted_pr[:limit_pagerank]):
                node_type = "Internal" if node.startswith(library_name) else "External"
                role = module_roles.get(node, "External Lib")
                short_name = node.replace(library_name + ".", "")
                lines.append(f"| {i+1} | `{short_name}` | {score:.4f} | {node_type} | {role} |")
            lines.append("\n")
        except Exception: pass

    lines.append("### 🗺️ Dependency & Architecture Map")
    mermaid_lines = ["graph TD"]
    mermaid_lines.append("    classDef core fill:#f96,stroke:#333,stroke-width:2px;")
    mermaid_lines.append("    classDef external fill:#9cf,stroke:#333,stroke-width:1px;")
    
    # YOU CAN MODIFY THIS TO LIMIT NODES IF NEEDED !!!!
    # Select top nodes to visualize: top 20
    # if HAS_NETWORKX: top_nodes = set(n for n, s in sorted_pr[:20])
    # else: top_nodes = set(x[0] for x in internal_modules_rank.most_common(20))
    # if not, we will show all
    # ⚠️[ARG6]   
    LIMIT_NODES = limit_dependency_graph
    if HAS_NETWORKX:
        top_nodes = set(n for n, s in sorted_pr[:LIMIT_NODES])
    else:
        top_nodes = set(x[0] for x in internal_modules_rank.most_common(LIMIT_NODES))
    
    # --- Use pure numeric ID mapping ---
    # 1. Collect all node names to be drawn
    nodes_to_map = set()
     
    # Collect nodes from dependency relationships
    source_data = G.edges() if HAS_NETWORKX else []
    if not HAS_NETWORKX:
        for src, targets in dependency_graph.items():
            for tgt in targets: source_data.append((src, tgt))
            
    dependency_edges = []
    for u, v in source_data:
        if u in top_nodes or v in top_nodes:
            # Simplify naming logic
            short_u = u.replace(library_name + ".", "").split('.')[-1]
            short_v = v.replace(library_name + ".", "").split('.')[-1]
            if not v.startswith(library_name): short_v = v.split('.')[0]
            
            if short_u == short_v: continue
            
            nodes_to_map.add(u)
            nodes_to_map.add(v)
            dependency_edges.append((u, v, short_u, short_v))

    # Collect nodes from inheritance relationships
    inheritance_edges = []
    all_classes = []
    for mod in submodules:
        for name, obj in inspect.getmembers(mod, inspect.isclass):
            if getattr(obj, "__module__", "").startswith(library_name):
                all_classes.append((name, obj))
                
    # YOU CAN MODIFY THIS TO LIMIT NODES IF NEEDED !!!! 
    # JUST MOVE THE DOUBLE FOR LOOP INSIDE THE IF CONDITION MARKDER TO LIMIT PROCESSING
    # ⚠️[ARG7]
    LIMIT_CLASSES = limit_inheritance_graph 
    if len(all_classes) < LIMIT_CLASSES:
        for name, obj in all_classes:
            for base in obj.__bases__:
                base_name = base.__name__ 
                if base_name == 'object': continue
                    
                nodes_to_map.add(name)
                nodes_to_map.add(base_name)
                
                base_module = base.__module__.split('.')[0]
                inheritance_edges.append((name, base_name, base_module))
    # 2. Build ID mapping table
    id_map = {name: f"id_{i}" for i, name in enumerate(nodes_to_map)}

    # 3. Draw dependency relationships
    edges_drawn = set()
    for u, v, label_u, label_v in dependency_edges:
        uid = id_map[u]
        vid = id_map[v]
        
        edge_key = f"{uid}->{vid}"
        if edge_key in edges_drawn: continue
        edges_drawn.add(edge_key)
        
        arrow = "-.->" if not v.startswith(library_name) else "-->"
        mermaid_lines.append(f'    {uid}["{label_u}"] {arrow} {vid}["{label_v}"]')
        
        if u.startswith(library_name): mermaid_lines.append(f"    class {uid} core;")
        else: mermaid_lines.append(f"    class {uid} external;")
        
        if v.startswith(library_name): mermaid_lines.append(f"    class {vid} core;")
        else: mermaid_lines.append(f"    class {vid} external;")

    # 4. Draw inheritance relationships
    for cls_name, base_name, base_mod in inheritance_edges:
        cid = id_map[cls_name]
        bid = id_map[base_name]
        
        mermaid_lines.append(f'    {cid}["{cls_name}"] ==> {bid}["{base_name}"]')
        mermaid_lines.append(f"    class {cid} core;")
        
        if base_mod != library_name:
            mermaid_lines.append(f"    class {bid} external;")
        else:
            mermaid_lines.append(f"    class {bid} core;")

    lines.append("```mermaid")
    lines.append("\n".join(mermaid_lines))
    lines.append("```\n")

    # ===============================================
    # Phase 1.5: Global Call Graph & Extraction Guide 
    # ===============================================
    lines.append("## 🚀 Global Execution Flow & Extraction Guide")
    lines.append("This graph visualizes how data flows between functions across the entire project.")
    
    global_call_graph, dependency_map = generate_global_call_graph(submodules, library_name)
    if global_call_graph:
        lines.append("```mermaid")
        lines.append(global_call_graph)
        lines.append("```\n")
        
        # --- Navigator: Extraction Guide ---
        lines.append("### ✂️ Navigator: Snippet Extractor")
        # ⚠️[ARG8]
        lines.append(f"Want to use a specific function without the whole library? Here is the **Dependency Closure** for **Top {limit_extraction_guide}** key functions.")
        
        # Select several important functions for analysis (if PageRank is available, select the top-ranked; otherwise, select the most called)
        # For simplicity, select the top 20 functions that appear most frequently as callers in the dependency_map
        top_funcs = sorted(dependency_map.keys(), key=lambda k: len(dependency_map[k]), reverse=True)[:limit_extraction_guide]
        
        if top_funcs:
            for func in top_funcs:
                closure = get_dependency_closure(func, dependency_map)
                lines.append(f"#### To extract `{func}`:")
                lines.append(f"> You need these **{len(closure)}** components:")
                lines.append(f"`{', '.join(sorted(list(closure)))}`")
                lines.append("")
        else:
            lines.append("_Not enough call data to generate extraction guide._")

    else:
        lines.append("_No internal function calls detected (or code structure is too dynamic)._\n")

    # ===============================================
    # Phase 2: Surface Level Inspection & Logic Flow
    # ===============================================
    lines.append("## 📑 Top-Level API Contents & Logic Flow")

    if hasattr(main_module, "__all__"):
        all_names = main_module.__all__
        using_all = True
    else:
        all_names = dir(main_module)
        using_all = False
    
    members_data = []

    for name in all_names:
        if not include_private and not using_all and name.startswith("_"):
            continue
        try: obj = getattr(main_module, name)
        except AttributeError: continue
        obj_module = getattr(obj, "__module__", None)
        is_imported = False
        if obj_module and not obj_module.startswith(library_name): is_imported = True
        if not include_imported and is_imported:
             if not using_all: continue
        members_data.append((name, obj, is_imported))

    classes = []
    functions = []
    for name, obj, is_imported in members_data:
        display_name = name + (" (imported)" if is_imported else "")
        if inspect.isclass(obj): classes.append((display_name, obj))
        elif inspect.isfunction(obj) or inspect.isbuiltin(obj): functions.append((display_name, obj))

    def get_info(obj):
        """
        Description
        ----------
        Extract signature and docstring summary for a Python object.

        Args
        -----
        obj : Any
            The function or class object to inspect.

        Returns
        --------
        Tuple[str, str]
            A tuple containing (signature_string, docstring_summary).

        Notes
        -------
        - Handles cases where signature extraction fails (e.g., built-ins) by falling back to `__text_signature__` or default.
        """
        try: sig = str(inspect.signature(obj))
        except (ValueError, TypeError):
            sig = getattr(obj, "__text_signature__", "(...)")
            if sig is None: sig = "(...)"
        doc = inspect.getdoc(obj) or "No documentation available."
        return sig, doc

    if functions:
        lines.append("### 🔧 Functions")
        for name, func in functions:
            sig, doc = get_info(func)
            lines.append(f"#### `{name}{sig}`")
            lines.append(f"> {doc.splitlines()[0] if doc else ''}")
            
            lines.append(f"<details><summary>Full Docstring</summary>\n\n```text\n{doc}\n```\n</details>\n")

            flow_chart = generate_function_flowchart(func)
            if flow_chart:
                lines.append("\n**Logic Flow:**")
                lines.append("```mermaid")
                lines.append(flow_chart)
                lines.append("```\n")

    if classes:
        lines.append("### 📦 Classes")
        for name, cls in classes:
            sig, doc = get_info(cls)
            lines.append(f"#### `class {name}{sig}`")
            lines.append(f"{doc.splitlines()[0] if doc else ''}\n")
            
            methods = inspect.getmembers(cls, predicate=lambda x: inspect.isfunction(x) or inspect.ismethod(x))
            if methods:
                lines.append("| Method | Signature | Description |")
                lines.append("| :--- | :--- | :--- |")
                for m_name, m_obj in methods:
                    if not include_private and m_name.startswith("_") and m_name != "__init__":
                        continue
                    m_sig, m_doc = get_info(m_obj)
                    short_doc = m_doc.splitlines()[0] if m_doc else "-"
                    short_doc = short_doc.replace("|", "\\|")
                    lines.append(f"| **{m_name}** | `{m_sig}` | {short_doc} |")
            lines.append("\n")

    # --- Output ---
    content = "\n".join(lines)
    
    if output_path:
        # 1. Save Markdown (original logic)
        md_path = output_path
        if not md_path.endswith(".md"):
            md_path += ".md"
        
        output_dir = os.path.dirname(md_path)
        if output_dir and not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir)
            except OSError as e:
                print(f"❌ Error creating directory {output_dir}: {e}")
                return

        try:
            with open(md_path, "w", encoding="utf-8") as f:
                f.write(content)
            print(f"✅ Markdown report saved to: {os.path.abspath(md_path)}")
            
            # 2. (New) Automatically generate HTML version
            html_path = md_path.replace(".md", ".html")
            html_content = convert_md_to_html(content, f"Analysis Report: {library_name}")
            
            with open(html_path, "w", encoding="utf-8") as f:
                f.write(html_content)
            print(f"📊 Interactive HTML report saved to: {os.path.abspath(html_path)}")
            print(f"   (Open the HTML file in your browser to see rendered charts)")
            
        except IOError as e:
            print(f"❌ Error writing file: {e}")
    else:
        print(content)