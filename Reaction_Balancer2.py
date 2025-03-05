import re
from sympy import Matrix, lcm

def chem_compound(compound):
    """Extract elements and their counts from the formulas."""
    pattern = r'([A-Z][a-z]*)(\d*)'  # Extracting elements and atom counts
    compound_match = re.findall(pattern, compound)  # Searching for matches
    compound_dict = {}  # Dictionary to store element and count
    for element, count in compound_match:
        count = int(count) if count else 1  # Convert count to int, default to 1
        compound_dict[element] = compound_dict.get(element, 0) + count
    return compound_dict

def chem_equation(equation):
    """Parse the equation into reactants and products."""
    reactants, products = equation.replace(' ', '').split('=')  # Remove spaces and split at '='
    reactants = reactants.split('+')  # Split reactants
    products = products.split('+')  # Split products
    return reactants, products  # Returns a tuple of two lists

def ver_bal(compounds_in_eq, bal_coeffs, reactants_count):
    """Verify that the equation is correctly balanced."""
    compound_dict = {}
    for i, compound in enumerate(compounds_in_eq):  # Loop through compounds
        parsed = chem_compound(compound)  # Extract element counts
        for element, count in parsed.items():
            if element not in compound_dict:
                compound_dict[element] = [0, 0]  # Initialize reactant and product counts
            if i < reactants_count:
                compound_dict[element][0] += bal_coeffs[i] * count  # Update reactant count
            else:
                compound_dict[element][1] += bal_coeffs[i] * count  # Update product count
    return all(reactant == product for reactant, product in compound_dict.values())

def balance_equation(equation):
    """Balance the given chemical equation."""
    reactants, products = chem_equation(equation)
    reactants_count = len(reactants)
    compounds_in_eq = reactants + products  # Merge reactants and products
    unique_elements = set()  # Collect unique elements
    compound_matrix = []  # Initialize matrix

    for compound in compounds_in_eq:
        parsed = chem_compound(compound)  # Extract element counts
        unique_elements.update(parsed.keys())  # Store unique elements

    unique_elements = sorted(unique_elements)  # Sort for consistency

    for element in unique_elements:
        row_matrix = [chem_compound(compound).get(element, 0) for compound in compounds_in_eq]
        compound_matrix.append(row_matrix)  # Append row to matrix

    X = Matrix(compound_matrix)  # Convert list into a mathematical matrix
    chem_coeffs = X.nullspace()

    if not chem_coeffs:  # Handle empty nullspace cases
        return "Error: Unable to balance"
    
    chem_coeffs = chem_coeffs[0]  # Extract the first solution vector

    # Ensure integer coefficients
    lcm_val = lcm([c.q if hasattr(c, 'q') else 1 for c in chem_coeffs])
    bal_coeffs = [abs(int(c * lcm_val)) for c in chem_coeffs]

    # Ensure coefficients are positive
    if any(c < 0 for c in bal_coeffs):
        bal_coeffs = [-c for c in bal_coeffs]

    # Verify if the equation is actually balanced
    if not ver_bal(compounds_in_eq, bal_coeffs, reactants_count):
        return "Error: Unable to balance"

    # Formatting the output
    reaction_part = ' + '.join(f"{bal_coeffs[i]}{reactants[i]}" for i in range(len(reactants)))
    prod_part = ' + '.join(f"{bal_coeffs[i + len(reactants)]}{products[i]}" for i in range(len(products)))

    return f"{reaction_part} = {prod_part}"

if __name__ == "__main__":
    equation = input("Enter a chemical equation (e.g., H2 + O2 = H2O): ")
    bal_eq = balance_equation(equation)
    print("Balanced Eq:", bal_eq)
