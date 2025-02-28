import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, DataStructs, rdMolDescriptors
import py3Dmol
import io

# 页面标题
st.title("Molecule Comparison")

# 创建左右两列用于对比两个分子
col_left, col_right = st.columns(2)


# 定义一个函数用于绘制 2D 和 3D 分子
def visualize_molecule(smiles_input, atom_indices_input, col):
    with col:
        st.subheader("Molecule Visualization")
        if smiles_input:
            try:
                mol = Chem.MolFromSmiles(smiles_input)
                if mol:
                    atom_indices = list(map(int, atom_indices_input.split(','))) if atom_indices_input else []
                    col_2d, col_3d = st.columns(2)

                    # **2D 可视化**
                    with col_2d:
                        st.subheader("2D Visualization")
                        Chem.AllChem.Compute2DCoords(mol)
                        img = Draw.MolToImage(mol, size=(300, 300), highlightAtoms=atom_indices)
                        img_buffer = io.BytesIO()
                        img.save(img_buffer, format="PNG")
                        img_buffer.seek(0)
                        st.image(img)

                    # **3D 可视化**
                    with col_3d:
                        st.subheader("3D Visualization")
                        mol_3d = Chem.AddHs(mol)
                        AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())
                        mol_block = Chem.MolToMolBlock(mol_3d)
                        view = py3Dmol.view(width=400, height=300)
                        view.addModel(mol_block, "mol")
                        view.setStyle({"stick": {}})
                        view.zoomTo()
                        for idx in atom_indices:
                            atom_pos = mol_3d.GetConformer().GetAtomPosition(idx)
                            view.addSphere({"center": {"x": atom_pos.x, "y": atom_pos.y, "z": atom_pos.z},
                                            "radius": 0.6, "color": "red", "opacity": 0.8})
                        view_html = view._make_html()
                        st.components.v1.html(view_html, height=300, width=400)
                else:
                    st.error("Invalid SMILES string. Please try again.")
            except Exception as e:
                st.error(f"Error processing SMILES: {e}")


# 左侧分子
with col_left:
    st.subheader("Molecule 1")
    smiles_input_1 = st.text_input("Enter SMILES String for Molecule 1", "CC(=O)Oc1ccccc1C(O)=O")
    atom_indices_input_1 = st.text_input("Enter list of atom indices to highlight for Molecule 1", "0,1,2")
    visualize_molecule(smiles_input_1, atom_indices_input_1, col_left)

# 右侧分子
with col_right:
    st.subheader("Molecule 2")
    smiles_input_2 = st.text_input("Enter SMILES String for Molecule 2", "C1=CC=CC(=C1C(O)=O)OC(C)=O")
    atom_indices_input_2 = st.text_input("Enter list of atom indices to highlight for Molecule 2", "0,1")
    visualize_molecule(smiles_input_2, atom_indices_input_2, col_right)

# 分子相似性和立体化学比较
st.subheader("Molecule Similarity and Stereochemistry Check")
if smiles_input_1 and smiles_input_2:
    try:
        mol1 = Chem.MolFromSmiles(smiles_input_1)
        mol2 = Chem.MolFromSmiles(smiles_input_2)

        if mol1 and mol2:
            smiles1_canonical = Chem.MolToSmiles(mol1, isomericSmiles=True)
            smiles2_canonical = Chem.MolToSmiles(mol2, isomericSmiles=True)
            are_equal_graph = Chem.MolToSmiles(mol1) == Chem.MolToSmiles(mol2)
            are_equal_stereo = smiles1_canonical == smiles2_canonical

            # 计算 Tanimoto 相似性
            fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
            fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
            tanimoto_sim = DataStructs.FingerprintSimilarity(fp1, fp2)

            # 计算基于 BRICS 的分子相似性
            fp1_brics = rdMolDescriptors.GetHashedBRICSFingerprint(mol1)
            fp2_brics = rdMolDescriptors.GetHashedBRICSFingerprint(mol2)
            brics_sim = DataStructs.DiceSimilarity(fp1_brics, fp2_brics)

            st.markdown("### Similarity Results")
            st.markdown(
                f"""
                | Comparison                        | Result              |
                |-----------------------------------|---------------------|
                | Molecular graphs identical?      | {'**Yes**' if are_equal_graph else '**No**'} |
                | Stereochemically identical?      | {'**Yes**' if are_equal_stereo else '**No**'} |
                | Tanimoto Similarity (0-1)       | **{tanimoto_sim:.4f}** |
                | BRICS Similarity (Dice) (0-1)   | **{brics_sim:.4f}** |
                """
            )

            st.markdown("### Canonical SMILES")
            st.markdown(
                f"""
                | Molecule   | Canonical SMILES                       |
                |------------|----------------------------------------|
                | Molecule 1 | `{smiles1_canonical}`                 |
                | Molecule 2 | `{smiles2_canonical}`                 |
                """
            )
        else:
            st.error("One or both SMILES strings are invalid.")
    except Exception as e:
        st.error(f"Error occurred while comparing molecules: {e}")
