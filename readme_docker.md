# Dockerise a shiny app on Google cloud

## Major steps

1. Make sure that the shiny application works on R studio
2. Build docker file 
  * from shiny-verse
  * Install packages as needed (load install_packages.R)
  * Use the shiny-server command
3. Test docker locally
4. Upload to GitHub, this will automatically force compilation on Google cloud and Docker site.

If the Docker application has not been yet uploaded to Google
1. Upload to Google cloud
  * Set up enough memory (2 Go for pr2 primers)
2. Upload to Docker web site

## Building docker image

* Do not use package renv which is really heavy....
* Start Powershell under windows

```
cd "C:/daniel.vaulot@gmail.com/Papers/2020 Egge MicroPolar 18S/data/micropolar-protists_datapaper"

# Buid with cache

docker build . -t micropolar-protists

docker run --rm -p 8080:8080 micropolar-protists

```

Test locally

* http://localhost:8080/


## Push to Cloud run


Utilize Google Builds to build image on the cloud

```
gcloud auth login

gcloud auth configure-docker

gcloud config set project tactile-bolt-247111

gcloud builds submit --tag asia.gcr.io/tactile-bolt-247111/micropolar-protists
```

Deploy to Google Cloud Run (Need only first time)

```
gcloud run deploy --image asia.gcr.io/tactile-bolt-247111/micropolar-protists --platform managed --max-instances 1
```

Effectuer ensuite un mappage de domaine sur:

http:/micropolar-protists.metapr2.org


## Push container to Docker repository

* Can also be done with Docker Desktop

```
docker images

docker tag micropolar-protists/micropolar-protists:v1.0.0

docker push vaulot/micropolar-protists:v1.0.0
```

## Docker misc

* List running containers

```
docker container ls
```

* Stop a container
```
docker stop 12a32e8928ef
```

* Remove dangling caches
```
docker builder prune
docker image prune
```

* Buid without cache
```
docker build --no-cache . -t asia.gcr.io/aerial-citron-246112/pr2-primers
```

* Management

```
# see All containers
docker ps -a

# See images and remove

docker images
docker rmi xxxxx

# Remove containers, image, cache and volumes

docker system prune --volumes
```

