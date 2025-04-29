require 'net/http'
require 'uri'
          
server="http://localhost:8000"
path = "/api/sequences/get_sequence_meta_data/PP760209"

# Construct the full URL with query parameters
url = URI.parse(server + path)


# Make the HTTP request
http = Net::HTTP.new(url.host, url.port)
request = Net::HTTP::Get.new(url, { 'Content-Type' => 'application/json' })

response = http.request(request)

puts response.body

if response.code != "200"
    puts "Invalid response: #{response.code}"
    puts response.body
    exit
end