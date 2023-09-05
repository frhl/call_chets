version="0.0.2"
echo "${version} [$(git rev-parse main | cut -c1-7)]" > .version
