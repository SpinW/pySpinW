import ast
import operator

""" Safe evaluation of generator expressions """

operators = {
    ast.Add: operator.add,
    ast.Sub: operator.sub,
    ast.Mult: operator.mul,
    ast.Div: operator.truediv,
    ast.USub: operator.neg,
}


def evaluate_syntax_tree_node(node):
    """ Recursive evaluation of a mathematical expression """

    if isinstance(node, ast.Expression):
        return evaluate_syntax_tree_node(node.body)

    elif isinstance(node, ast.BinOp):
        left = evaluate_syntax_tree_node(node.left)
        right = evaluate_syntax_tree_node(node.right)
        return operators[type(node.op)](left, right)

    elif isinstance(node, ast.UnaryOp):
        operand = evaluate_syntax_tree_node(node.operand)
        return operators[type(node.op)](operand)

    elif isinstance(node, ast.Constant):  # For Python 3.8+
        return node.value

    else:
        raise ValueError(f"Unsupported operation: {node}")


def evaluate_algebra(string_expr: str):
    """ Evaluate a algebraic expression """

    try:
        parsed_expr = ast.parse(string_expr, mode='eval')
        return evaluate_syntax_tree_node(parsed_expr.body)

    except Exception as e:
        raise ValueError(f"Invalid expression: {string_expr}") from e

