import re
import numpy as np #useful for matrix operations
from sympy import Matrix, lcm #lcm ensures there are integer coefficients

def chem_compound(compound): #using this to extract elements and their counts from the formulas
    pattern = r'([A-Z][a-z]*)(\d*)' #extracting the elements and their atom counts
    compound_match = re.findall(pattern, compound) # searching the compound and matching them to the pattern
    compound_dict = {} # an empty dictionary that will store the element and its count
    for element, count in compound_match:
        count = int(count) if count else 1 #changes the count in the tuple to an integer, if theres no number it changes it to 1
        compound_dict[element] = compound_dict.get(element, 0) + count #see molecular weight calculator files for explanation
    return compound_dict

#parsing the equation
def chem_equation(equation): #parsing the equation into its reactants and products
    reactants, products = equation.replace(' ', '').split('=') #repalces all spaces with no spaces, splits the equation at the equal to sign
    reactants = reactants.split('+') #splits reactants at each plus sign so we have individual components
    products = products.split('+') #same as line above
    return reactants, products #the defined fuction is a tuple of two lists, reactants and products

def ver_bal(compounds_in_eq, bal_coeffs, reactants_count): #making sure the equation balances, and doesn't spit out the nullspace solution
    compound_dict = {}

    for i, compound in enumerate(compounds_in_eq): #looping through all compounds
        parsed = chem_compound(compound) #gets the element counts
        for element, count in parsed.items():
            if element not in compound_dict:
                compound_dict[element] = [0,0] #setting reactant and product count to 0 if element is not present
            if i < reactants_count:
                compound_dict[element][0] += bal_coeffs[i] * count #if compound is a reactant, update count in first part of the list 
            else:
                compound_dict[element][1] += bal_coeffs[i] * count #if compound is a product, update count in second part of the list
    return all(reactant == product for reactant, product in compound_dict.values())
            

#Creating the matrix to balance the equation
def balance_equation(equation):
    reactants, products = chem_equation(equation)
    reactants_count = len(reactants)
    compounds_in_eq = reactants + products #merges the lists 
    unique_elements = set() #collects the unique elements in the equation
    compound_matrix = [] #empty list
    for compound in compounds_in_eq:
        parsed = chem_compound(compound) #for each compound in the equation, we are using them to populate the dictionary
        unique_elements.update(parsed.keys()) #updating the unique elements set with the keys from the dictionary
    unique_elements = sorted(unique_elements) #returns a list

    for element in unique_elements:
        row_matrix = [chem_compound(compound).get(element, 0) for compound in compounds_in_eq]
        #in the above line, for each compound found in the merged compounds_in_eq list, we are using chem_compound(compound) to create a dictionary for that compound
        #we are then then looking for each element in "elements" from the compound's dictionary and returning the value
        #process is repeated for every compound in the equation
        
        compound_matrix.append(row_matrix) #compound_matrix becomes a list of lists where each row is for one element, each column is for one compound
    X = Matrix(compound_matrix) #converts the list into a mathematical matrix
    chem_coeffs = X.nullspace() #we are looking for a vector when multiplied by the matrix X, the result is zero vector.
    if not chem_coeffs:
        return "Error: Unable to balance"
    chem_coeffs = chem_coeffs[0] #[0] means we are taking the first solution from the list of column vectors which are possible solutions
 
    
    lcm_val = lcm([c.q if hasattr(c, 'q') else 1 for c in chem_coeffs]) #gets the lcm of the denominators
    bal_coeffs = [abs(int(c*lcm_val)) for c in chem_coeffs] #multiplying the lcm by the coeffs, returns a list of balanced coeffs

    if any(c < 0 for c in bal_coeffs): #ensuring coefficients are integers and positive
        bal_coeffs = [-c for c in bal_coeffs]

    #this checks if the equation is actually balanced
    if not ver_bal(compounds_in_eq, bal_coeffs, reactants_count):
        return "Error: Unable to balance"

    #formatting the output
    reaction_part = '+' .join(f"{bal_coeffs[i]}{reactants[i]}" for i in range (len(reactants))) #i in length of reactants,
    #if there is 2 reactants, i is 1 or 2
    prod_part = '+' .join(f"{bal_coeffs[i+len(reactants)]}{products[i]}" for i in range (len(products)))
    return f"{reaction_part} = {prod_part}"

if __name__ == "__main__":
    equation = input("Enter a chemical equation (eg. H2 + O2 = H2O):")
    bal_eq = balance_equation(equation)
    print ("Balanced Eq:", bal_eq)

        


    
    
    
