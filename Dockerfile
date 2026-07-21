FROM rust:1.86-alpine3.21 AS builder

WORKDIR /app

RUN apk add --no-cache musl-dev libc-dev

COPY Cargo.toml Cargo.lock ./
RUN cargo fetch --locked

COPY . .

RUN cargo build --release --locked

FROM alpine:3.21

WORKDIR /app
COPY --from=builder /app/target/release/krakenclip /app/krakenclip

ENTRYPOINT ["/app/krakenclip"]
CMD ["--help"] 