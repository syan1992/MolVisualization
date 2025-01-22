import streamlit as st
import pandas as pd
import os
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource
from bokeh.models import HoverTool
from st_aggrid import AgGrid, GridOptionsBuilder

dataset_folder = "Datasets"  # 替换为你的数据集所在文件夹路径

# 获取数据集名称列表
dataset_files = [
    'BACE', 'BBBP', 'CLINTOX', 'ESOL', 'Freesolv', 'HIV', 'Lipophilicity', 'Sider', 'Tox21', 'PIN1'
]
# 页面布局
col1, col2 = st.columns([1, 3])  # 左1/右3比例分割

# 左侧：显示数据集名称列表
with col1:
    selected_dataset = st.radio(
        "Select a dataset:",
        dataset_files,
        index=0
    )

# 右侧：显示选中数据集的统计结果
with col2:
    st.header(selected_dataset)
    dataset_path = os.path.join(dataset_folder, selected_dataset + '.csv')

    # 加载数据集
    try:
        data = pd.read_csv(dataset_path)
        # 创建选项卡导航栏
        tab1, tab2, tab3 = st.tabs(["Preview", "Statistics", "Visualization"])

        # 选项卡 1：数据预览
        with tab1:
            st.write("数据预览：")
            num_rows_to_display = st.slider("选择要显示的行数：", min_value=5, max_value=50, value=10, step=5)
            preview_data = data.head(num_rows_to_display)

            # 配置 AgGrid 样式和交互
            gb = GridOptionsBuilder.from_dataframe(preview_data)
            gb.configure_default_column(cellStyle={'textAlign': 'center'})  # 设置列居中
            grid_options = gb.build()

            AgGrid(preview_data, gridOptions=grid_options, height=300, fit_columns_on_grid_load=True)
        # 选项卡 2：数据统计
        with tab2:
            st.write("数据统计：")
            stats_data = data.describe()
            stats_data = stats_data.reset_index()  # 将索引变为普通列
            stats_data.rename(columns={'index': 'metrics'}, inplace=True)
            # 配置 AgGrid 的样式和功能
            gb = GridOptionsBuilder.from_dataframe(stats_data)
            gb.configure_default_column(
                cellStyle={'textAlign': 'center'}  # 单元格内容居中
            )
            gb.configure_grid_options(domLayout='normal')  # 默认布局
            grid_options = gb.build()

            # 使用 AgGrid 显示数据统计表
            AgGrid(
                stats_data,
                gridOptions=grid_options,
                height=300,
                fit_columns_on_grid_load=True  # 列宽自适应
            )
        with tab3:
            st.write("分类数据可视化：")

            # 定义统计列规则
            classification_columns = {
                "BACE": ["Class"],
                "BBBP": ["p_np"],
                "CLINTOX": data.columns[1:3].tolist(),
                "Sider": data.columns[1:28].tolist(),
                "Tox21": data.columns[:12].tolist(),
                "HIV": ["HIV_active"],
            }

            regression_columns = {
                "Lipophilicity": ['exp'],
                "ESOL": ['ESOL predicted log solubility in mols per litre'],
                "Freesolv": ['expt'],
                "PIN1": ['Potency', 'Efficacy', 'P_E']
            }

            if selected_dataset in classification_columns:
                columns_to_plot = classification_columns[selected_dataset]

                for i, column in enumerate(columns_to_plot):
                    # 如果是第一列或每五个图后，创建一个新的列布局
                    if i % 5 == 0:
                        cols = st.columns(5)

                    # 获取当前列的值分布
                    class_counts = data[column].value_counts()
                    source = ColumnDataSource(data={
                        'Class': class_counts.index.astype(str),
                        'Counts': class_counts.values
                    })

                    # 使用 Bokeh 创建柱状图
                    hover = HoverTool(tooltips=[
                        ("Class", "@Class"),
                        ("Counts", "@Counts")
                    ])

                    p = figure(
                        x_range=list(class_counts.index.astype(str)),  # 转换为标准列表
                        plot_height=300,
                        plot_width=250,
                        title=f"分类数据分布: {column}",
                        toolbar_location=None,
                        tools=[hover]
                    )
                    p.vbar(
                        x='Class', top='Counts', width=0.5, source=source, color="#90EE90"  # 浅蓝色
                    )
                    p.xgrid.grid_line_color = None
                    p.y_range.start = 0
                    p.xaxis.axis_label = "Class"
                    p.yaxis.axis_label = "Counts"

                    # 将图表放置到当前列中
                    cols[i % 5].bokeh_chart(p)

            elif selected_dataset in regression_columns:
                columns_to_plot = regression_columns[selected_dataset]

                for column in columns_to_plot:
                    # 使用 Bokeh 绘制分布图
                    p = figure(
                        plot_height=400,
                        plot_width=600,
                        title=f"回归数据分布: {column}",
                        toolbar_location=None,
                        tools=""
                    )

                    # 计算直方图数据
                    hist, edges = pd.cut(data[column], bins=40, retbins=True)
                    counts = hist.value_counts().sort_index()

                    # 修正：将 counts 和 edges 转换为列表并添加到数据源
                    counts_list = list(counts)  # 将 pd.Series 转换为列表
                    left_edges = list(edges[:-1])  # 左边界
                    right_edges = list(edges[1:])  # 右边界

                    # 准备数据源
                    source = ColumnDataSource(data={
                        'left': left_edges,  # 左边界
                        'right': right_edges,  # 右边界
                        'top': counts_list  # 对应的高度
                    })

                    # 绘制分布图
                    p.quad(
                        bottom=0,
                        top='top',
                        left='left',
                        right='right',
                        source=source,
                        fill_color="#90EE90",
                        line_color="white"
                    )
                    p.xaxis.axis_label = column
                    p.yaxis.axis_label = "Counts"

                    st.bokeh_chart(p)
            else:
                st.warning("当前数据集没有分类或回归列可统计。")
    except Exception as e:
        st.error(f"加载数据集时出错: {e}")
