version="0.0.3"
echo "${version} [$(git rev-parse main | cut -c1-7)]" > .version
