<template>
  <v-container class="text-center" id="home">
    <img id="profilepic" src="../assets/profilepic.jpg" />
    <br />
    <h1>{{ greeting }}!</h1>
    <AboutMe />
  </v-container>
</template>

<script>
import AboutMe from "@/components/AboutMe.vue";
export default {
  components: {
    AboutMe
  },
  data: function() {
    return {
      greeting: "Hello world, it's me"
    };
  }
};
</script>

<style scoped>
#home {
  margin: auto;
  padding: 5%;
}
h1 {
  padding: 2%;
  font-size: 2em;
  text-align: center;
}
#profilepic {
  border-radius: 50%;
  width: 300px;
  padding: 10px;
}
</style>
