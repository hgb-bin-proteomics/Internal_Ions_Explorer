# Dockerfile for Internal Ion Explorer
# author: Micha Birklbauer
# version: 1.3.0

FROM python:3.12

LABEL maintainer="micha.birklbauer@gmail.com"

RUN mkdir internal_ions
COPY ./ internal_ions/
WORKDIR internal_ions

RUN pip install --upgrade pip
RUN pip install --upgrade setuptools
RUN pip install --no-cache-dir uv

RUN uv sync
RUN uv cache clean

CMD  ["sh", "-c", "uv run -- streamlit run streamlit_app.py"]
