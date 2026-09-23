"""
Decompose a reaction's gene-protein-reaction (GPR) rule into its candidate
enzyme complexes -- the minimal gene sets that can each, independently,
satisfy the rule.

`map_gpr` (in `mteapy.utils`) collapses a whole GPR straight to a single
expression score (OR -> max/absmax, AND -> min). This module stops one step
earlier: it expands the GPR to disjunctive-normal-form (an OR of ANDs) and
returns each AND-clause as one complex, so a reaction with an OR in its GPR
can be scored per-complex rather than as one flattened number. This is the
basis for context-aware task scoring: two tissues that both "run" the same
reaction may do so via two different complexes (an isozyme swap), which a
single flattened GPR score cannot distinguish.
"""
import re

import boolean


class UnsupportedGPRExpression(Exception):
    def __init__(self, expression):
        self.message = f"Could not parse GPR expression into complexes: {expression!r}"
        super().__init__(self.message)


def _get_enzymes_recursive(expression, or_operator: str = "|", and_operator: str = "&"):
    list_of_enzymes = []
    if expression.isliteral:
        list_of_enzymes.append(str(expression))
    else:
        expression = expression.simplify()
        if expression.operator == or_operator:
            for arg in expression.args:
                if arg.isliteral:
                    list_of_enzymes.append([str(arg)])
                else:
                    list_of_enzymes.extend(_get_enzymes_recursive(arg))
        elif expression.operator == and_operator:
            enzymes = []
            simple = True
            for arg in expression.args:
                subunits = _get_enzymes_recursive(arg)
                if len(subunits) > 1:
                    simple = False
                    enzymes = [enzymes + s for s in subunits]
                else:
                    enzymes.append(subunits[0])
            if simple:
                list_of_enzymes.append(enzymes)
            else:
                list_of_enzymes.extend(enzymes)
        else:
            raise UnsupportedGPRExpression(expression)

    return list_of_enzymes


def get_enzymes(gene_reaction_rule: str) -> tuple[tuple[str, ...], ...]:
    """
    Expand a GPR rule to disjunctive-normal-form and return each AND-clause
    as one candidate enzyme complex -- a minimal set of genes that can, on
    its own, satisfy the rule. Complexes are deduplicated and each complex's
    genes are sorted, so e.g. `"A or (A and B)"` correctly collapses to the
    single minimal complex `("A",)`, not `("A",), ("A", "B")`.

    Parameters
    ----------
    gene_reaction_rule: str
        A reaction's GPR rule as a plain boolean-logic string, e.g.
        `model.reactions.get_by_id(rxn_id).gene_reaction_rule`.

    Returns
    -------
    complexes: tuple[tuple[str, ...], ...]
        One tuple of gene IDs per candidate complex. A reaction with no GPR
        (`""`) or a single gene returns one complex with one gene.
    """
    if not gene_reaction_rule:
        return ()

    # Gene IDs containing "-" would otherwise be parsed as boolean NOT by
    # the `boolean` library; swap to "_" for parsing, then swap back.
    dash_found = bool(re.search("-", gene_reaction_rule))
    rule = re.sub("-", "_", gene_reaction_rule) if dash_found else gene_reaction_rule

    algebra = boolean.BooleanAlgebra()
    expression = algebra.parse(rule, simplify=True)

    if expression.isliteral:
        list_of_enzymes = [[str(expression)]]
    else:
        expression = expression.simplify()
        list_of_enzymes = _get_enzymes_recursive(expression)
        wrong_parse = any(
            not isinstance(subunit, str) for enzyme in list_of_enzymes for subunit in enzyme
        )
        if wrong_parse:
            # The fast path above assumes the GPR is already close to OR-of-
            # ANDs form; a genuinely nested mix (e.g. "A and (B or C) and D")
            # needs full distributive expansion to reach true DNF.
            expression = expression.distributive()
            list_of_enzymes = _get_enzymes_recursive(expression)

    if dash_found:
        for enzyme in list_of_enzymes:
            for i, gene in enumerate(enzyme):
                enzyme[i] = re.sub("_", "-", gene)

    return tuple(sorted(set(tuple(sorted(enzyme)) for enzyme in list_of_enzymes)))
