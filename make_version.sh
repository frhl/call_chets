version="0.1.8"
echo "${version} [$(git rev-parse main | cut -c1-7)]" > .version
