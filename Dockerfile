# Dockerfile for Internal Ion Explorer
# author: Micha Birklbauer
# version: 1.3.0

FROM python:3.14

LABEL maintainer="micha.birklbauer@gmail.com"

RUN mkdir app
COPY ./ app/
WORKDIR app

RUN pip install --upgrade pip
RUN pip install --upgrade setuptools
RUN pip install --no-cache-dir uv

RUN uv sync --no-dev --no-cache

CMD  ["sh", "-c", "uv run --no-sync streamlit run streamlit_app.py"]
