#!/bin/bash
# Streamlit Dashboard 启动脚本

echo "================================"
echo "  ChronoDecon Dashboard"
echo "================================"
echo ""

# 切换到仓库根目录（脚本所在目录）
REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$REPO_ROOT"

# 检查是否在仓库根目录
if [ ! -f "app.py" ] || [ ! -d "chrono_decon" ]; then
    echo "错误: 请在仓库根目录下运行此脚本（app.py 和 chrono_decon/ 应同目录）"
    echo "当前目录: $(pwd)"
    exit 1
fi

# 使用 miniconda3 的 Python 环境（有 matchms）
PYTHON_ENV="/home/knan/miniconda3/bin/python3"
STREAMLIT_CMD="/home/knan/miniconda3/bin/streamlit"

echo "Python 环境: $PYTHON_ENV"
echo "Streamlit 命令: $STREAMLIT_CMD"
echo ""
echo "验证 matchms 安装..."
$PYTHON_ENV -c "import matchms; print(f'✅ matchms version: {matchms.__version__}')" 2>/dev/null || {
    echo "❌ matchms 未安装或不可访问"
    echo "尝试安装 matchms..."
    $PYTHON_ENV -m pip install matchms || {
        echo "❌ 无法安装 matchms"
        exit 1
    }
}

echo ""
echo "启动 Streamlit Dashboard..."
echo "访问地址: http://localhost:8501"
echo ""
echo "按 Ctrl+C 停止服务器"
echo ""

cd "$REPO_ROOT"
$STREAMLIT_CMD run app.py
