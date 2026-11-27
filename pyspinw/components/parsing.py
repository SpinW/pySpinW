""" Parsing of polarisation description"""

import ast
import operator

from .components import (
    Sperp,
    Sxx, Sxy, Sxz,
    Syx, Syy, Syz,
    Szx, Szy, Szz,
    Mab,
    Mxx, Mxy, Mxz,
    Myx, Myy, Myz,
    Mzx, Mzy, Mzz,
)


BIN_OPS = {
    ast.Add: operator.add,
    ast.Sub: operator.sub,
    ast.Mult: operator.mul,
    ast.Div: operator.truediv,
}

UNARY_OPS = {
    ast.UAdd: operator.pos,
    ast.USub: operator.neg,
}


def evaluate(node, variables):
    "Parse and algebraic expression "
    if isinstance(node, ast.Expression):
        return evaluate(node.body, variables)

    if isinstance(node, ast.Constant):
        if isinstance(node.value, (int, float, complex)):
            return node.value
        raise ValueError("Invalid constant")

    if isinstance(node, ast.Name):
        if node.id not in variables:
            raise ValueError(f"Unknown variable: {node.id}")
        return variables[node.id]

    if isinstance(node, ast.BinOp):
        op = BIN_OPS.get(type(node.op))
        if op is None:
            raise ValueError("Unsupported operator")

        left = evaluate(node.left, variables)
        right = evaluate(node.right, variables)
        return op(left, right)

    if isinstance(node, ast.UnaryOp):
        op = UNARY_OPS.get(type(node.op))
        if op is None:
            raise ValueError("Unsupported unary operator")

        return op(evaluate(node.operand, variables))

    raise ValueError(f"Unsupported expression: {type(node).__name__}")


def parse_polarisation(expression):


    variables = (
        Sperp,
        Sxx, Sxy, Sxz,
        Syx, Syy, Syz,
        Szx, Szy, Szz,
        Mab,
        Mxx, Mxy, Mxz,
        Myx, Myy, Myz,
        Mzx, Mzy, Mzz,
    )

    tree = ast.parse(expression, mode="eval")
    return evaluate(tree, variables)