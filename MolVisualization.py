import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import py3Dmol
import io

# 页面标题
#st.title("Molecule Visualization")

# SMILES 输入框
smiles_input = st.text_input("Enter SMILES String", "CCO")  # 默认输入乙醇 (ethanol)

# 原子序号输入框
atom_indices_input = st.text_input("Enter list of atom indices to highlight (comma-separated)", "0,1,2")  # 默认输入0,1,2

if smiles_input:
    try:
        # 将 SMILES 转换为 RDKit 分子对象
        mol = Chem.MolFromSmiles(smiles_input)

        # 检查分子是否有效
        if mol:
            # 将输入的原子序号转换为整数列表
            atom_indices = list(map(int, atom_indices_input.split(',')))

            # 创建两个列
            col1, col2 = st.columns(2)

            # **2D 可视化**
            with col1:
                st.subheader("2D Visualization")
                # 为 2D 图高亮原子
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

                # 提供下载按钮
                st.download_button(
                    label="Download 2D Image",
                    data=img_buffer,
                    file_name="molecule_2d.png",
                    mime="image/png",
                )

            # **3D 可视化**
            with col2:
                st.subheader("3D Visualization")
                # 使用 RDKit 生成 3D 坐标
                mol_3d = Chem.AddHs(mol)  # 添加氢原子
                AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())  # 生成 3D 坐标

                # 将分子对象转换为 MOL 文件字符串
                mol_block = Chem.MolToMolBlock(mol_3d)

                # 使用 py3Dmol 显示 3D 结构
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
                view_html = view._make_html()  # 确保 3D 内容生成 HTML
                st.components.v1.html(view_html, height=300, width=400)
        else:
            st.error("Invalid SMILES string. Please try again.")
    except Exception as e:
        st.error(f"Error occurred while processing the SMILES string: {e}")
