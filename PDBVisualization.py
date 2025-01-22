import streamlit as st
import py3Dmol

# 页面标题
st.title("PDB Visualization")

# 上传 PDB 文件
uploaded_file = st.file_uploader("Upload a PDB file", type=["pdb"])

# 如果用户未上传文件，使用默认文件
if not uploaded_file:
    st.info("No file uploaded. Using default PDB file: Datasets/1pin.pdb")
    pdb_name = "1pin.pdb"  # 默认文件名称
    try:
        with open("Datasets/1pin.pdb", "r") as f:
            pdb_content = f.read()
    except FileNotFoundError:
        st.error("Default file not found. Please upload a PDB file.")
        pdb_content = None
else:
    # 读取用户上传的文件
    pdb_name = uploaded_file.name  # 上传文件名称
    pdb_content = uploaded_file.read().decode("utf-8")

# 可视化样式选项
visualization_styles = ["cartoon", "stick", "sphere", "surface"]
selected_style = st.selectbox("Select Visualization Style", visualization_styles, index=0)

# 如果有 PDB 内容，进行可视化
if pdb_content:
    try:
        # 使用 py3Dmol 可视化 PDB
        view = py3Dmol.view(width=800, height=600)
        view.addModelsAsFrames(pdb_content)  # 加载 PDB 数据
        view.setStyle({'model': -1}, {selected_style: {'color': 'lightgreen'}})  # 根据选择设置样式
        view.zoomTo()  # 自动调整视角

        # 渲染视图
        view_html = view._make_html()
        st.components.v1.html(view_html, height=600, width=800)

        # 显示 PDB 文件名称
        st.subheader(f"PDB:{pdb_name.split('.')[0]}")
    except Exception as e:
        st.error(f"Error occurred while processing the PDB file: {e}")
