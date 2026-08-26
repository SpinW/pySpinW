""" The contents of this file is used for parsing text forms of symmetry operations, e.g. (y-x,x,z+1/2)"""

import ast
import operator

OPS = {
    ast.Add: operator.add,
    ast.Sub: operator.sub,
    ast.Div: operator.truediv,
    ast.Pow: operator.pow,
    ast.USub: operator.neg,
}

def evaluate(expression, variables):
    """ Numerical evaluation of algebraic expression written as a string

    :param expression: String form expression, e.g. x+y
    :param variables: dictionary of string/value pairs e.g. {"x": 1, "y": 2}

    """
    tree = ast.parse(expression, mode="eval")

    def eval_node(node):
        if isinstance(node, ast.Constant):
            return node.value

        if isinstance(node, ast.Name):
            return variables[node.id]

        if isinstance(node, ast.BinOp):
            op = OPS[type(node.op)]
            return op(eval_node(node.left), eval_node(node.right))

        if isinstance(node, ast.UnaryOp):
            op = OPS[type(node.op)]
            return op(eval_node(node.operand))

        raise ValueError(f"Unsupported expression: {ast.dump(node)}")

    return eval_node(tree.body)

if __name__ == "__main__":
    print(evaluate("x+y-z",{"x": 1, "y": 3, "z": 10}))
