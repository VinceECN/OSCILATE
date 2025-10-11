SymPy basics
============

SymPy is the Python symbolic calculus library on which the package relies. 
Below are some basics about the possibilities it offers. 
For more information, refer to `SymPy's doc <https://docs.sympy.org/latest/index.html>`_.

Symbols and functions
---------------------

To declare a parameter :math:`t` and a function :math:`f(t)` using SymPy, run::

    from sympy import symbols, Function
    t = symbols("t")
    f = Function("f")(t)

Working with expressions
------------------------

Combinations of symbols and/or Functions form SymPy expressions. For instance, ::

    expr = 2 f

is an expression. 
Expressions can be simplified, developed or factored using the ``.simplify()``, ``.expand()``, ``.factor()`` methods.

One can also collect specific terms in a (developed) expression using ::

    dic = expr.collect(term, evaluate=False)

This returns a dictionary whose keys are powers of ``term`` and values are the factor of these keys in ``expr``.

Display
-------

The lines ::

    from sympy.physics.vector.printing import init_vprinting
    init_vprinting(use_latex=True, forecolor='White') # Initialise latex printing 

allow to initialise the output display such that

- LaTeX is used for rendering if the python environment supports it.

- The LaTeX colour is set to white (ideal for dark modes)

Additionally, one can display the unformatted LaTeX expression (for copy-paste in a LaTeX document for instance) using ::

    from sympy.physics.vector.printing import vlatex
    print(vlatex(expr))

