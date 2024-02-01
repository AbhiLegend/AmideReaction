import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from io import BytesIO

# Function to convert SMILES to molecule and render image
def render_molecule(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            img = Draw.MolToImage(mol)
            byte_io = BytesIO()
            img.save(byte_io, format='PNG')
            return byte_io.getvalue(), mol
        else:
            return None, None
    except Exception as e:
        st.error(f"Error in rendering molecule: {e}")
        return None, None

# Amide coupling function
def amide_coupling(smile1, smile2):
    mol1 = Chem.MolFromSmiles(smile1)
    mol2 = Chem.MolFromSmiles(smile2)
    smarts = "[C:1](=[O:2])O.[Nh:3] >> [C:1](=[O:2])[Nh:3]"
    rxn3 = AllChem.ReactionFromSmarts(smarts)
    products = rxn3.RunReactants([mol1, mol2])
    resulting_smile_list = []
    try:
        for i in range(len(products)):
            resulting_smile = Chem.MolToSmiles(products[i][0])
            resulting_smile_list.append(resulting_smile)
    except:
        pass
    return resulting_smile_list

# Streamlit app
def main():
    st.title("Amide Coupling Reaction Visualization")

    # User input for SMILES strings
    smile1 = st.text_input("Enter SMILES string for molecule 1 (e.g., a carboxylic acid)")
    smile2 = st.text_input("Enter SMILES string for molecule 2 (e.g., an amine)")

    if smile1 and smile2:
        # Render reactants
        image_data1, mol1 = render_molecule(smile1)
        image_data2, mol2 = render_molecule(smile2)

        if image_data1 and image_data2:
            col1, col2 = st.columns(2)
            with col1:
                st.image(image_data1, caption="Molecule 1", use_column_width=True)
            with col2:
                st.image(image_data2, caption="Molecule 2", use_column_width=True)

            # Perform reaction and visualize products
            product_smiles = amide_coupling(smile1, smile2)
            if product_smiles:
                st.write("Products:")
                for smi in product_smiles:
                    prod_img, _ = render_molecule(smi)
                    st.image(prod_img, caption=f"Product: {smi}", use_column_width=True)
            else:
                st.error("No products generated. Please check the SMILES strings.")
        else:
            st.error("Unable to render one or both molecules. Please check the SMILES strings.")

if __name__ == "__main__":
    main()
