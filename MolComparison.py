import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
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
                # 将 SMILES 转换为 RDKit 分子对象
                mol = Chem.MolFromSmiles(smiles_input)

                # 检查分子是否有效
                if mol:
                    # 将输入的原子序号转换为整数列表
                    if atom_indices_input:
                        atom_indices = list(map(int, atom_indices_input.split(',')))
                    else:
                        atom_indices = []

                    # 创建上下两部分
                    col_2d, col_3d = st.columns(2)

                    # **2D 可视化**
                    with col_2d:
                        st.subheader("2D Visualization")
                        Chem.AllChem.Compute2DCoords(mol)  # 生成 2D 坐标

                        if len(atom_indices) > 0:
                            img = Draw.MolToImage(mol, size=(300, 300), highlightAtoms=atom_indices)  # 高亮指定原子
                        else:
                            img = Draw.MolToImage(mol, size=(300, 300))

                        # 将图像保存到 BytesIO
                        img_buffer = io.BytesIO()
                        img.save(img_buffer, format="PNG")
                        img_buffer.seek(0)

                        st.image(img)

                    # **3D 可视化**
                    with col_3d:
                        st.subheader("3D Visualization")
                        mol_3d = Chem.AddHs(mol)  # 添加氢原子
                        AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())  # 生成 3D 坐标

                        mol_block = Chem.MolToMolBlock(mol_3d)

                        view = py3Dmol.view(width=400, height=300)
                        view.addModel(mol_block, "mol")  # 加载 MOL 文件字符串
                        view.setStyle({"stick": {}})  # 使用 Stick 样式
                        view.zoomTo()  # 自动调整视图

                        # 高亮 3D 图中的原子
                        for idx in atom_indices:
                            atom_pos = mol_3d.GetConformer().GetAtomPosition(idx)  # 获取原子坐标
                            view.addSphere(
                                {"center": {"x": atom_pos.x, "y": atom_pos.y, "z": atom_pos.z},
                                 "radius": 0.6, "color": "red", "opacity": 0.8}  # 高亮原子
                            )

                        # 渲染 3D 视图
                        view_html = view._make_html()
                        st.components.v1.html(view_html, height=300, width=400)
                else:
                    st.error("Invalid SMILES string. Please try again.")
            except Exception as e:
                st.error(f"Error occurred while processing the SMILES string: {e}")

# 左侧分子
with col_left:
    st.subheader("Molecule 1")
    smiles_input_1 = st.text_input("Enter SMILES String for Molecule 1", "CC(=O)Oc1ccccc1C(O)=O")  # 默认输入乙醇
    atom_indices_input_1 = st.text_input("Enter list of atom indices to highlight for Molecule 1", "0,1,2")
    visualize_molecule(smiles_input_1, atom_indices_input_1, col_left)

# 右侧分子
with col_right:
    st.subheader("Molecule 2")
    smiles_input_2 = st.text_input("Enter SMILES String for Molecule 2", "C1=CC=CC(=C1C(O)=O)OC(C)=O")  # 默认输入环己烷
    atom_indices_input_2 = st.text_input("Enter list of atom indices to highlight for Molecule 2", "0,1")
    visualize_molecule(smiles_input_2, atom_indices_input_2, col_right)

# 分子相似性和立体化学比较
st.subheader("Molecule Similarity and Stereochemistry Check")
if smiles_input_1 and smiles_input_2:
    try:
        mol1 = Chem.MolFromSmiles(smiles_input_1)
        mol2 = Chem.MolFromSmiles(smiles_input_2)

        if mol1 and mol2:
            # 标准化 SMILES 表示
            smiles1_canonical = Chem.MolToSmiles(mol1, isomericSmiles=True)
            smiles2_canonical = Chem.MolToSmiles(mol2, isomericSmiles=True)

            # 比较分子图是否相同
            are_equal_graph = Chem.MolToSmiles(mol1) == Chem.MolToSmiles(mol2)

            # 比较立体化学是否相同
            are_equal_stereo = smiles1_canonical == smiles2_canonical

            # 使用 Markdown 优化显示
            st.markdown("### Similarity Results")
            st.markdown(
                f"""
                | Comparison                        | Result              |
                |-----------------------------------|---------------------|
                | Molecular graphs identical?      | {'**Yes**' if are_equal_graph else '**No**'} |
                | Stereochemically identical?      | {'**Yes**' if are_equal_stereo else '**No**'} |
                """
            )

            # 显示规范化的 SMILES
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
