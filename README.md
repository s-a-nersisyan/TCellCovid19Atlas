# TCellCovid19Atlas
## DB
- docker-compose up -d (run docker image)
- docker-compose stop (stop container, but do not remove it)
- docker-compose down (stop container and remove it)
- docker ps (see running containers)
- docker exec -it containername /bin/bash (attach to container's shell)
- psql -h localhost -d TCellCovid19 -U covid19admin (connect to db)
